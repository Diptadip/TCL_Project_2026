clear all
close all
clc

global TCL

%% ── 0.  Load mechanistic model parameters ────────────────────────────────
load TCL_MechModel_Parameters.mat   

A  = TCL.phy;       
B  = TCL.gama;      
C  = TCL.C_mat;     
D  = TCL.D_mat;     
Ts = TCL.Samp_T;    

[n_op, n_st] = size(C);
[~,    n_ip] = size(B);

Ys = TCL.Ys(:);   
Us = TCL.Us(:);   

fprintf('=== TCL State-Tracking MPC (Mechanistic) ===\n');
fprintf('  State dim n_st = %d,  Inputs n_ip = %d,  Outputs n_op = %d\n', n_st, n_ip, n_op);
fprintf('  Operating point:  T1_ss = %.2f degC   T2_ss = %.2f degC\n', Ys(1), Ys(2));
fprintf('  Operating point:  H1_ss = %.2f %%      H2_ss = %.2f %%\n\n', Us(1), Us(2));

%% ── 1.  MPC Tuning Parameters ────────────────────────────────────────────
Np    = 15;                           
Nc    = 5;                            

% Weights (Penalizing states)
Wx    = diag([15,  1 ]);             % state tracking weight (n_st x n_st)
Wdelu = diag([1.0, 1.0]);             % input-move weight     (n_ip x n_ip)

fprintf('=== MPC Tuning ===\n');
fprintf('  Np = %d,  Nc = %d,  Ts = %g s\n', Np, Nc, Ts);
fprintf('  Wx    = diag([%.1f %.1f])\n', Wx(1,1), Wx(2,2));
fprintf('  Wdelu = diag([%.1f %.1f])\n\n', Wdelu(1,1), Wdelu(2,2));

%% ── 2.  Kalman Filter Design ─────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st);
R_kf = TCL.R;                        
[~, ~, L_kf] = dare(A', C', Q_kf, R_kf);
L_kf = L_kf';                        

fprintf('Kalman gain L_kf computed (DARE).\n\n');

%% ── 3.  Steady-State Target Computation ─────────────────────────────────
M_ss     = [(eye(n_st) - A), -B ; C, zeros(n_op, n_ip)];
M_ss_inv = pinv(M_ss);

fprintf('Steady-state target gain M_ss_inv pre-computed.\n\n');

%% ── 4.  Pre-compute QP Matrices (State Prediction) ──────────────────────
A_aug = [A, B; zeros(n_ip, n_st), eye(n_ip)];
B_aug = [B; eye(n_ip)];
% State extractor matrix
Cx_aug = [eye(n_st), zeros(n_st, n_ip)]; 

n_aug = n_st + n_ip;

Psi_x   = zeros(Np*n_st, n_aug);
Theta_x = zeros(Np*n_st, Nc*n_ip);

A_aug_pow = eye(n_aug);
for i = 1:Np
    A_aug_pow         = A_aug_pow * A_aug;
    rows_i            = (i-1)*n_st + 1 : i*n_st;
    Psi_x(rows_i, :)  = Cx_aug * A_aug_pow;

    for j = 1:min(i, Nc)
        rows_ij  = rows_i;
        cols_j   = (j-1)*n_ip + 1 : j*n_ip;
        A_aug_ij = eye(n_aug);
        for ii = 1:i-j
            A_aug_ij = A_aug_ij * A_aug;
        end
        Theta_x(rows_ij, cols_j) = Cx_aug * A_aug_ij * B_aug;
    end
end

Wx_bar    = kron(eye(Np), Wx);
Wdelu_bar = kron(eye(Nc), Wdelu);

% Hessian Matrix (H) for quadprog
H_qp = Theta_x' * Wx_bar * Theta_x + Wdelu_bar;
H_qp = (H_qp + H_qp') / 2; % Ensure exact symmetry

% Constraint Matrix M 
M_mat  = kron(tril(ones(Nc)), eye(n_ip));
A_cons = [M_mat; -M_mat]; 

qp_options = optimoptions('quadprog', 'Display', 'off');

%% ── 5.  Experiment Timing ───────────────────────────────────────────────
N        = 550;
k_warmup = 30;
k_settle = 100;
k_step1  = 101;
k_step2  = 301;

delta_r_step1 = [ +4; -2 ];
delta_r_step2 = [ -2; +4 ];

kT    = (0:N-1)';
t_sec = kT * Ts;

u_H = [100; 100] - Us;    
u_L = [  0;   0] - Us;    

%% ── 6.  Filter Parameters ───────────────────────────────────────────────
beta_r = 0.95;   phy_r = beta_r * eye(n_op);
alfa_e = 0.95;   phy_e = alfa_e * eye(n_op);

%% ── 7.  Plant Simulation Choice ─────────────────────────────────────────
Plant_Sim = 1;

if Plant_Sim == 0
    fprintf('Plant: NONLINEAR  (ode45 / TCL_Dynamics)\n\n');
else
    fprintf('Plant: LINEAR mechanistic  (TCL.phy, TCL.gama)\n\n');
end

%% ── 8.  Run simulation for Noise_ON = 0 and 1 ───────────────────────────
for Noise_ON = 0:1

    if Noise_ON
        fprintf('\n====  Noisy Simulation  ====\n\n');
        rng(42);
        vk = mvnrnd(zeros(n_op,1), TCL.R, N)';
        fig_offset = 4;
        run_label  = 'Noisy';
    else
        fprintf('\n====  Noise-Free Simulation  ====\n\n');
        vk = zeros(n_op, N);
        fig_offset = 0;
        run_label  = 'Noise-Free';
    end

    % ── Allocate ─────────────────────────────────────────────────────────
    xk_hat = zeros(n_st, N);     
    uk     = zeros(n_ip, N);     
    yk     = zeros(n_op, N);     
    rk_f   = zeros(n_op, N);     
    ek_f   = zeros(n_op, N);     
    xs_k   = zeros(n_st, N);     
    us_k   = zeros(n_ip, N);     
    xk_true= zeros(n_st, N);     

    u_prev = zeros(n_ip, 1);

    % ── Main loop ─────────────────────────────────────────────────────────
    for k = 1:N

        if k < k_step1, rk = zeros(n_op, 1);
        elseif k < k_step2, rk = delta_r_step1;
        else, rk = delta_r_step2; end

        if k == 1, rk_f(:,k) = (eye(n_op) - phy_r) * rk;
        else, rk_f(:,k) = phy_r * rk_f(:,k-1) + (eye(n_op) - phy_r) * rk; end

        raw_ek = rk_f(:,k) - yk(:,k);
        if k == 1, ek_f(:,k) = (eye(n_op) - phy_e) * raw_ek;
        else, ek_f(:,k) = phy_e * ek_f(:,k-1) + (eye(n_op) - phy_e) * raw_ek; end

        % Steady-state targets
        rhs       = [zeros(n_st,1); rk_f(:,k)];
        xus       = M_ss_inv * rhs;
        xs_k(:,k) = xus(1:n_st);
        us_k(:,k) = xus(n_st+1:end);

        if k <= k_warmup
            uk(:,k)     = zeros(n_ip,1);
            u_prev      = zeros(n_ip,1);
            xk_hat(:,k) = zeros(n_st,1);
        else
            x_aug = [xk_hat(:,k); u_prev];

            % STATE Tracking Reference (x_s)
            X_ref = repmat(xs_k(:,k), Np, 1);   

            % Linear term for QP: f = Theta' * W * (Psi*x - Ref)
            f_qp = Theta_x' * Wx_bar * (Psi_x * x_aug - X_ref);

            % Update Inequality Bounds (dynamic based on u_prev)
            b_cons = [repmat(u_H - u_prev, Nc, 1); 
                      repmat(-u_L + u_prev, Nc, 1)];

            % Solve Constrained QP
            [du_full, ~, exitflag] = quadprog(H_qp, f_qp, A_cons, b_cons, [], [], [], [], [], qp_options);

            % Robustness fallback
            if isempty(du_full) || exitflag < 0
                du_opt = zeros(n_ip, 1);
            else
                du_opt = du_full(1:n_ip);
            end

            uk(:,k) = u_prev + du_opt;
            u_prev  = uk(:,k);
        end

        if k < N
            if Plant_Sim == 0
                TCL.Uk = TCL.Us + uk(:,k);
                [~, Xt] = ode45('TCL_Dynamics', [0 Ts], TCL.Xs + xk_true(:,k));
                xk_true(:,k+1) = Xt(end,:)' - TCL.Xs;
            else
                xk_true(:,k+1) = A * xk_true(:,k) + B * uk(:,k);
            end
            yk(:,k+1) = C * xk_true(:,k+1) + vk(:,k+1);

            y_pred        = C * xk_hat(:,k);
            inno          = yk(:,k) - y_pred;
            xk_hat(:,k+1) = A * xk_hat(:,k) + B * uk(:,k) + L_kf * inno;
        end

    end 

    % ── Absolute arrays ───────────────────────────────────────────────────
    Yk   = yk   + repmat(Ys, 1, N);
    Uk   = uk   + repmat(Us, 1, N);
    Rk_f = rk_f + repmat(Ys, 1, N);
    Xs_k = xs_k + repmat(TCL.Xs, 1, N);
    Us_k = us_k + repmat(Us,    1, N);

    Uk_ss1 = mean(Uk(:, k_step2-20 : k_step2-1), 2);
    Uk_ss2 = mean(Uk(:, N-19 : N),               2);
    fprintf('Simulated SS heater values:\n');
    fprintf('  Step1: H1=%.1f%%  H2=%.1f%%\n', Uk_ss1(1), Uk_ss1(2));
    fprintf('  Step2: H1=%.1f%%  H2=%.1f%%\n', Uk_ss2(1), Uk_ss2(2));

    % ── FIGURE A ────────────────────────────────────────
    figure(fig_offset+1)
    clf
    subplot(211)
    plot(kT, Yk(1,:), 'b-', 'LineWidth',1.5), hold on
    plot(kT, Rk_f(1,:), 'b--', 'LineWidth',1.0)
    grid on, ylabel('T_1  (deg C)')
    title(['State-Tracking MPC – ' run_label ' – Output Profiles'])
    legend('Y_1(k)', 'R_1(k)', 'Location','Best')
    xline(k_step1-1, 'k:', 'Step 1', 'LabelVerticalAlignment','bottom')
    xline(k_step2-1, 'k:', 'Step 2', 'LabelVerticalAlignment','bottom')

    subplot(212)
    plot(kT, Yk(2,:), 'r-', 'LineWidth',1.5), hold on
    plot(kT, Rk_f(2,:), 'r--', 'LineWidth',1.0)
    grid on, xlabel('Sample k'), ylabel('T_2  (deg C)')
    legend('Y_2(k)', 'R_2(k)', 'Location','Best')
    xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

    % ── FIGURE B ─────────────────────────────────────────
    figure(fig_offset+2)
    clf
    subplot(211)
    stairs(kT, Uk(1,:)', 'b-', 'LineWidth',1.5), hold on
    yline(5,  'k--', 'LineWidth',0.8), yline(80, 'k--', 'LineWidth',0.8)
    grid on, ylabel('U_1(k)  (% heater)')
    title(['State-Tracking MPC – ' run_label ' – Heater Inputs'])
    legend('U_1(k)', 'Safety bounds', 'Location','Best')
    xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

    subplot(212)
    stairs(kT, Uk(2,:)', 'r-', 'LineWidth',1.5), hold on
    yline(5,  'k--', 'LineWidth',0.8), yline(80, 'k--', 'LineWidth',0.8)
    grid on, xlabel('Sample k'), ylabel('U_2(k)  (% heater)')
    legend('U_2(k)', 'Safety bounds', 'Location','Best')
    xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

    % ── FIGURE C ─────────────────────────────────────────
    figure(fig_offset+3)
    clf
    subplot(211)
    stairs(kT, ek_f(1,:)', 'b-', 'LineWidth',1.5), hold on
    yline(0, 'k--', 'LineWidth',0.8), grid on
    ylabel('e_{f,1}(k)  (deg C)')
    title(['State-Tracking MPC – ' run_label ' – Filtered Error e_f(k)'])
    xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

    subplot(212)
    stairs(kT, ek_f(2,:)', 'r-', 'LineWidth',1.5), hold on
    yline(0, 'k--', 'LineWidth',0.8), grid on
    xlabel('Sample k'), ylabel('e_{f,2}(k)  (deg C)')
    xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

    % ── FIGURE D ────────────────────────────
    figure(fig_offset+4)
    clf
    subplot(211)
    plot(kT, Xs_k(1,:), 'b-', 'LineWidth',1.5), hold on
    plot(kT, Xs_k(2,:), 'r-', 'LineWidth',1.5)
    grid on, ylabel('x_s(k)  (deg C)')
    title(['State-Tracking MPC – ' run_label ' – SS Target States x_s(k)'])
    legend('x_{s,1}(k)', 'x_{s,2}(k)', 'Location','Best')
    xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

    subplot(212)
    stairs(kT, Us_k(1,:)', 'b-', 'LineWidth',1.5), hold on
    stairs(kT, Us_k(2,:)', 'r-', 'LineWidth',1.5)
    grid on, xlabel('Sample k'), ylabel('u_s(k)  (% heater)')
    title(['State-Tracking MPC – ' run_label ' – SS Target Inputs u_s(k)'])
    legend('u_{s,1}(k)', 'u_{s,2}(k)', 'Location','Best')
    xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

    if Noise_ON == 0
        save TCL_MPC_Mech_SimResults_NoiseFree  kT t_sec Yk Uk Rk_f ek_f Xs_k Us_k
        fprintf('Saved: TCL_MPC_Mech_SimResults_NoiseFree.mat\n');
    else
        save TCL_MPC_Mech_SimResults_Noisy      kT t_sec Yk Uk Rk_f ek_f Xs_k Us_k
        fprintf('Saved: TCL_MPC_Mech_SimResults_Noisy.mat\n');
    end
end 

fprintf('\nAll simulations complete.\n');