% =========================================================================
%  TCL_MPC_Mech_Simulation.m
%
%  PURPOSE
%  -------
%  Design and simulate a MIMO STATE-SPACE MPC controller for the
%  Temperature Control Lab (TCL) using the MECHANISTIC linearised
%  discrete-time state-space model (TCL.phy, TCL.gama).
%
%  MPC SCHEME  –  State-Space (SS-MPC)
%  ------------------------------------
%  The QP cost penalises STATE deviations from a steady-state target xs
%  and input increments Du:
%
%    J = sum_{i=1}^{Np} (x(k+i) - xs)' Wx (x(k+i) - xs)
%      + sum_{j=0}^{Nc-1} Du(k+j)' Wdelu Du(k+j)
%
%  Prediction model (perturbation form, delta-u):
%    x_aug(k+1) = A_aug * x_aug(k) + B_aug * Du(k)
%  where  x_aug = [x_hat ; u_prev]
%
%  Reference vector over horizon built from SS target state xs_k:
%    Xs_vec = [xs_k; xs_k; ... (Np times)]   (using state rows of Psi only)
%
%  Key difference vs Output-Tracking MPC:
%    SS-MPC  : reference = xs_k  (state target, from SS equations)
%    OT-MPC  : reference = rk_f  (output setpoint, repeated over horizon)
%
%  CONTROLLER STRUCTURE
%  --------------------
%    - Augmented state  x_aug(k) = [x_hat(k); u(k-1)]
%    - Psi_x   : maps x_aug to predicted STATE trajectory  (Np*n_st x n_aug)
%    - Theta_x : maps Delta_U to predicted STATE trajectory (Np*n_st x Nc*n_ip)
%    - Unconstrained MPC gain K_mpc_ss computed from state cost
%    - Kalman filter for state estimation (DARE, same R, Q as OT file)
%    - SS targets (xs_k, us_k) solved online from filtered setpoint rk_f
%    - Anti-windup input clamping to [0, 100]%
%
%  MPC TUNING PARAMETERS
%  ---------------------
%    Prediction horizon  Np    = 15  samples  (60 s lookahead)
%    Control horizon     Nc    = 5   samples  (20 s free moves)
%    State weight        Wx    = diag([15, 15])   (state error penalty)
%    Input-move weight   Wdelu = diag([1.0, 1.0]) (move suppression)
%    Kalman Q_kf         = 1e-3 * eye(n_st)
%    Kalman R_kf         = TCL.R   (from PRBS experiment)
%
%  FILTERS  (same as all project servo files)
%    Setpoint filter:    beta_r = 0.95
%    Innovation filter:  alfa_e = 0.95
%
%  EXPERIMENT TIMING
%    Ts       = 4 s,   N = 550 samples  (~36.7 min)
%    k_warmup = 30,    k_settle = 100
%    k_step1  = 101,   k_step2  = 301
%    Step 1 : delta_r = [+6; -4] deg C
%    Step 2 : delta_r = [-5; +7] deg C
%
%  PLANT SIMULATION
%    Plant_Sim = 0  ->  ode45 nonlinear ODE  (TCL_Dynamics.m)
%    Plant_Sim = 1  ->  linear mechanistic model  (TCL.phy, TCL.gama)
%
%  SIMULATIONS
%    Noise_ON = 0  ->  noise-free
%    Noise_ON = 1  ->  noisy (measurement noise from TCL.R)
%    Both run sequentially.
%
%  PLOTS (per simulation run)
%    Fig 1/5 – Yi(k) vs k and Ri(k) vs k
%    Fig 2/6 – Ui(k) vs k  (heater inputs, stairs)
%    Fig 3/7 – ef(k) vs k  (filtered error)
%    Fig 4/8 – xs(k) vs k  and  us(k) vs k  (SS target trajectories)
%
%  DEPENDENCIES
%    TCL_MechModel_Parameters.mat  (TCL struct)
%    TCL_Dynamics.m                (nonlinear ODE, needed if Plant_Sim=0)
%
%  OUTPUT
%    TCL_MPC_Mech_SimResults_NoiseFree.mat
%    TCL_MPC_Mech_SimResults_Noisy.mat
%      Variables (suffix _mech): Yk_mech, Uk_mech, Rk_f_mech,
%                                ek_f_mech, Xs_k_mech, Us_k_mech
% =========================================================================

clear all
close all
clc

global TCL

%% ── 0.  Load mechanistic model parameters ────────────────────────────────
load TCL_MechModel_Parameters.mat   % TCL struct

% Discrete-time mechanistic model matrices (perturbation form)
A  = TCL.phy;       % Ad  (n_st x n_st)
B  = TCL.gama;      % Bd  (n_st x n_ip)
C  = TCL.C_mat;     % Cd  (n_op x n_st)   = eye(2) for TCL
D  = TCL.D_mat;     % Dd  (n_op x n_ip)   = zeros(2)
Ts = TCL.Samp_T;    % 4 s

[n_op, n_st] = size(C);
[~,    n_ip] = size(B);

% Operating point
Ys = TCL.Ys(:);    % 2x1  steady-state temperatures [deg C]
Us = TCL.Us(:);    % 2x1  steady-state heater inputs [%]
Xs = TCL.Xs(:);    % 2x1  steady-state states        [deg C]

fprintf('=== TCL STATE-SPACE MPC (Mechanistic Model) ===\n');
fprintf('  MPC scheme : SS-MPC  (state error cost, xs reference)\n');
fprintf('  n_st=%d  n_ip=%d  n_op=%d   Ts=%g s\n', n_st, n_ip, n_op, Ts);
fprintf('  Op. point:  T1_ss=%.2f degC   T2_ss=%.2f degC\n', Ys(1), Ys(2));
fprintf('  Op. point:  H1_ss=%.2f %%     H2_ss=%.2f %%\n\n', Us(1), Us(2));

%% ── 1.  MPC Tuning Parameters ────────────────────────────────────────────
Np    = 15;                       % prediction horizon (samples)
Nc    = 5;                        % control horizon   (samples)
Wx    = diag([15,  15 ]);         % STATE error weight   (n_st x n_st)
Wdelu = diag([1.0, 1.0]);         % input-move weight    (n_ip x n_ip)

fprintf('=== MPC Tuning ===\n');
fprintf('  Np=%d  Nc=%d  Ts=%g s\n', Np, Nc, Ts);
fprintf('  Wx (state weight)    = diag([%.1f %.1f])\n', Wx(1,1), Wx(2,2));
fprintf('  Wdelu (move weight)  = diag([%.1f %.1f])\n\n', Wdelu(1,1), Wdelu(2,2));

%% ── 2.  Kalman Filter Design ─────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st);
R_kf = TCL.R;                    % measurement noise covariance (from PRBS)

[~, ~, L_kf] = dare(A', C', Q_kf, R_kf);
L_kf = L_kf';                   % n_st x n_op  observer gain

fprintf('Kalman gain L_kf computed (DARE).\n\n');

%% ── 3.  Steady-State Target Computation ─────────────────────────────────
%  Solve for xs, us given output setpoint perturbation delta_r:
%    (I-A)*xs - B*us = 0
%          C*xs      = delta_r
%  => M_ss * [xs; us] = [0; delta_r]

M_ss     = [(eye(n_st) - A), -B ; C, zeros(n_op, n_ip)];
M_ss_inv = pinv(M_ss);

fprintf('Steady-state target solver M_ss_inv pre-computed.\n\n');

%% ── 4.  Augmented State-Space (delta-u formulation) ─────────────────────
%  x_aug(k+1) = A_aug * x_aug(k) + B_aug * Du(k)
%  where x_aug = [x_hat(k); u(k-1)],  Du = u(k) - u(k-1)

n_aug = n_st + n_ip;
A_aug = [A,                 B             ;
         zeros(n_ip, n_st), eye(n_ip)     ];
B_aug = [B        ;
         eye(n_ip)];

% NOTE: For SS-MPC we build Psi_x (predicts STATES, not outputs)
%       C_aug maps x_aug -> STATE portion only:  C_x_aug = [eye(n_st), 0]
C_x_aug = [eye(n_st), zeros(n_st, n_ip)];   % n_st x n_aug

%% ── 5.  Build Prediction Matrices (STATE-based) ──────────────────────────
%  Predicted state trajectory:
%    X_pred = Psi_x * x_aug(k) + Theta_x * Delta_U
%  where X_pred = [x(k+1); x(k+2); ...; x(k+Np)]   size Np*n_st x 1

Psi_x   = zeros(Np * n_st,  n_aug);
Theta_x = zeros(Np * n_st,  Nc * n_ip);

A_aug_pow = eye(n_aug);
for i = 1 : Np
    A_aug_pow       = A_aug_pow * A_aug;
    rows_i          = (i-1)*n_st + 1 : i*n_st;
    Psi_x(rows_i,:) = C_x_aug * A_aug_pow;

    for j = 1 : min(i, Nc)
        cols_j = (j-1)*n_ip + 1 : j*n_ip;
        A_ij   = eye(n_aug);
        for ii = 1 : i-j
            A_ij = A_ij * A_aug;
        end
        Theta_x(rows_i, cols_j) = C_x_aug * A_ij * B_aug;
    end
end

%% ── 6.  Pre-compute SS-MPC Gain K_mpc_ss ────────────────────────────────
%  State cost over horizon:
%    J_x  = (X_pred - Xs_vec)' Wx_bar (X_pred - Xs_vec)
%         = DU' Theta_x' Wx_bar Theta_x DU
%           + 2*(Psi_x*x_aug - Xs_vec)' Wx_bar Theta_x DU + const
%
%  Unconstrained optimum:
%    DU* = (Theta_x' Wx_bar Theta_x + Wdelu_bar)^{-1}
%           * Theta_x' Wx_bar * (Xs_vec - Psi_x * x_aug)
%
%  Only the FIRST block row DU*(1) is applied (receding horizon).

Wx_bar    = kron(eye(Np), Wx);        % Np*n_st x Np*n_st
Wdelu_bar = kron(eye(Nc), Wdelu);     % Nc*n_ip x Nc*n_ip

H_ss      = Theta_x' * Wx_bar * Theta_x + Wdelu_bar;
H_ss      = (H_ss + H_ss') / 2;              % enforce symmetry
K_mpc_ss  = H_ss \ (Theta_x' * Wx_bar);      % Nc*n_ip x Np*n_st
K_mpc_ss  = K_mpc_ss(1:n_ip, :);             % first move block [n_ip x Np*n_st]

fprintf('SS-MPC gain K_mpc_ss computed  (size %dx%d).\n', ...
        size(K_mpc_ss,1), size(K_mpc_ss,2));
fprintf('  Reference vector : xs_k repeated Np times over state horizon.\n\n');

%% ── 7.  Experiment Timing & Setpoint Schedule ────────────────────────────
N        = 550;
k_warmup = 30;
k_settle = 100;
k_step1  = 101;
k_step2  = 301;

delta_r_step1 = [ +6; -4];
delta_r_step2 = [ -5; +7];

kT    = (0 : N-1)';
t_sec = kT * Ts;

u_H = [100; 100] - Us;    % upper perturbation bound
u_L = [  0;   0] - Us;    % lower perturbation bound

fprintf('=== Experiment Timing ===\n');
fprintf('  N=%d samples,  Total=%.1f min\n', N, N*Ts/60);
fprintf('  Warm-up k=1..%d,  Settle k=%d..%d\n', k_warmup, k_warmup+1, k_settle);
fprintf('  Step1   k=%d..%d,  Step2 k=%d..%d\n\n', k_step1, k_step2-1, k_step2, N);

%% ── 8.  Filter Parameters ────────────────────────────────────────────────
beta_r = 0.95;   phy_r = beta_r * eye(n_op);
alfa_e = 0.95;   phy_e = alfa_e * eye(n_op);

%% ── 9.  Plant Simulation Choice ─────────────────────────────────────────
%  Plant_Sim = 0  ->  ode45 nonlinear ODE  (TCL_Dynamics.m)
%  Plant_Sim = 1  ->  linear mechanistic discrete model
Plant_Sim = 1;

if Plant_Sim == 0
    fprintf('Plant: NONLINEAR  (ode45 / TCL_Dynamics)\n\n');
else
    fprintf('Plant: LINEAR mechanistic  (TCL.phy, TCL.gama)\n\n');
end

%% ── 10.  Run simulation for Noise_ON = 0 and 1 ──────────────────────────
for Noise_ON = 0:1

    if Noise_ON
        fprintf('\n====  Noisy Simulation  ====\n\n');
        rng(42);
        vk         = mvnrnd(zeros(n_op,1), TCL.R, N)';
        fig_offset = 4;
        run_label  = 'Noisy';
    else
        fprintf('\n====  Noise-Free Simulation  ====\n\n');
        vk         = zeros(n_op, N);
        fig_offset = 0;
        run_label  = 'Noise-Free';
    end

    % ── Allocate ─────────────────────────────────────────────────────────
    xk_hat  = zeros(n_st, N);
    xk_true = zeros(n_st, N);
    uk      = zeros(n_ip, N);
    yk      = zeros(n_op, N);
    rk_f    = zeros(n_op, N);
    ek_f    = zeros(n_op, N);
    xs_k    = zeros(n_st, N);
    us_k    = zeros(n_ip, N);

    u_prev  = zeros(n_ip, 1);

    % ── Main closed-loop loop ─────────────────────────────────────────────
    for k = 1 : N

        % ---- Setpoint schedule ------------------------------------------
        if k < k_step1
            rk = zeros(n_op, 1);
        elseif k < k_step2
            rk = delta_r_step1;
        else
            rk = delta_r_step2;
        end

        % ---- Setpoint filter --------------------------------------------
        if k == 1
            rk_f(:,k) = (eye(n_op) - phy_r) * rk;
        else
            rk_f(:,k) = phy_r * rk_f(:,k-1) + (eye(n_op) - phy_r) * rk;
        end

        % ---- Error and innovation filter --------------------------------
        raw_ek = rk_f(:,k) - yk(:,k);
        if k == 1
            ek_f(:,k) = (eye(n_op) - phy_e) * raw_ek;
        else
            ek_f(:,k) = phy_e * ek_f(:,k-1) + (eye(n_op) - phy_e) * raw_ek;
        end

        % ---- Compute SS targets from filtered setpoint ------------------
        rhs       = [zeros(n_st, 1); rk_f(:,k)];
        xus       = M_ss_inv * rhs;
        xs_k(:,k) = xus(1 : n_st);       % target state  (perturbation)
        us_k(:,k) = xus(n_st+1 : end);   % target input  (perturbation)

        % ---- Warm-up: heaters OFF ----------------------------------------
        if k <= k_warmup
            uk(:,k)     = zeros(n_ip, 1);
            u_prev      = zeros(n_ip, 1);
            xk_hat(:,k) = zeros(n_st, 1);

        else
            % ---- Build augmented state ----------------------------------
            x_aug = [xk_hat(:,k); u_prev];

            % ---- SS-MPC: reference = xs_k repeated Np times over STATE horizon
            %  Xs_vec is [xs_k; xs_k; ...(Np)] in STATE space (Np*n_st x 1)
            Xs_vec = repmat(xs_k(:,k), Np, 1);

            % ---- Optimal first move (SS-MPC law) ------------------------
            %  Du* = K_mpc_ss * (Xs_vec - Psi_x * x_aug)
            du_opt  = K_mpc_ss * (Xs_vec - Psi_x * x_aug);
            u_new   = u_prev + du_opt;
            u_new   = min(max(u_new, u_L), u_H);   % anti-windup clamp
            uk(:,k) = u_new;
            u_prev  = u_new;
        end

        % ---- Plant update -----------------------------------------------
        if k < N
            if Plant_Sim == 0
                TCL.Uk = Us + uk(:,k);
                [~, Xt] = ode45('TCL_Dynamics', [0, Ts], Xs + xk_true(:,k));
                xk_true(:,k+1) = Xt(end,:)' - Xs;
            else
                xk_true(:,k+1) = A * xk_true(:,k) + B * uk(:,k);
            end
            yk(:,k+1) = C * xk_true(:,k+1) + vk(:,k+1);

            % ---- Kalman filter update -----------------------------------
            y_pred        = C * xk_hat(:,k);
            inno          = yk(:,k) - y_pred;
            xk_hat(:,k+1) = A * xk_hat(:,k) + B * uk(:,k) + L_kf * inno;
        end

    end  % k loop

    % ── Build absolute arrays ─────────────────────────────────────────────
    Yk_mech   = yk   + repmat(Ys, 1, N);
    Uk_mech   = uk   + repmat(Us, 1, N);
    Rk_f_mech = rk_f + repmat(Ys, 1, N);
    ek_f_mech = ek_f;
    Xs_k_mech = xs_k + repmat(Xs, 1, N);
    Us_k_mech = us_k + repmat(Us, 1, N);

    % ── Performance metrics ───────────────────────────────────────────────
    ISE_T1 = sum((Yk_mech(1,:) - Rk_f_mech(1,:)).^2) * Ts;
    ISE_T2 = sum((Yk_mech(2,:) - Rk_f_mech(2,:)).^2) * Ts;
    IAE_T1 = sum(abs(Yk_mech(1,:) - Rk_f_mech(1,:))) * Ts;
    IAE_T2 = sum(abs(Yk_mech(2,:) - Rk_f_mech(2,:))) * Ts;

    fprintf('=== Tracking Performance [%s] ===\n', run_label);
    fprintf('  ISE T1 = %8.2f   |   ISE T2 = %8.2f\n', ISE_T1, ISE_T2);
    fprintf('  IAE T1 = %8.2f   |   IAE T2 = %8.2f\n', IAE_T1, IAE_T2);

    % ── FIGURE 1 – Output profiles ────────────────────────────────────────
    figure(fig_offset + 1)
    clf
    subplot(2,1,1)
    plot(kT, Yk_mech(1,:), 'b-', 'LineWidth', 1.5), hold on
    plot(kT, Rk_f_mech(1,:), 'b--', 'LineWidth', 1.0)
    grid on
    ylabel('T_1  (deg C)')
    title(['SS-MPC Mech – ' run_label ' – Output Profiles'])
    legend('Y_1(k)', 'R_1(k)', 'Location', 'Best')
    xline(k_step1-1, 'k:', 'Step 1', 'LabelVerticalAlignment', 'bottom')
    xline(k_step2-1, 'k:', 'Step 2', 'LabelVerticalAlignment', 'bottom')

    subplot(2,1,2)
    plot(kT, Yk_mech(2,:), 'r-', 'LineWidth', 1.5), hold on
    plot(kT, Rk_f_mech(2,:), 'r--', 'LineWidth', 1.0)
    grid on
    xlabel('Sample k')
    ylabel('T_2  (deg C)')
    legend('Y_2(k)', 'R_2(k)', 'Location', 'Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── FIGURE 2 – Input profiles ─────────────────────────────────────────
    figure(fig_offset + 2)
    clf
    subplot(2,1,1)
    stairs(kT, Uk_mech(1,:)', 'b-', 'LineWidth', 1.5), hold on
    yline(5,  'k--', 'LineWidth', 0.8)
    yline(80, 'k--', 'LineWidth', 0.8)
    grid on
    ylabel('U_1(k)  (% heater)')
    title(['SS-MPC Mech – ' run_label ' – Heater Inputs'])
    legend('U_1(k)', 'Safety bounds', 'Location', 'Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(2,1,2)
    stairs(kT, Uk_mech(2,:)', 'r-', 'LineWidth', 1.5), hold on
    yline(5,  'k--', 'LineWidth', 0.8)
    yline(80, 'k--', 'LineWidth', 0.8)
    grid on
    xlabel('Sample k')
    ylabel('U_2(k)  (% heater)')
    legend('U_2(k)', 'Safety bounds', 'Location', 'Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── FIGURE 3 – Filtered error ─────────────────────────────────────────
    figure(fig_offset + 3)
    clf
    subplot(2,1,1)
    stairs(kT, ek_f_mech(1,:)', 'b-', 'LineWidth', 1.5), hold on
    yline(0, 'k--', 'LineWidth', 0.8), grid on
    ylabel('e_{f,1}(k)  (deg C)')
    title(['SS-MPC Mech – ' run_label ' – Filtered Error e_f(k)'])
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(2,1,2)
    stairs(kT, ek_f_mech(2,:)', 'r-', 'LineWidth', 1.5), hold on
    yline(0, 'k--', 'LineWidth', 0.8), grid on
    xlabel('Sample k')
    ylabel('e_{f,2}(k)  (deg C)')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── FIGURE 4 – SS targets xs(k) and us(k) ────────────────────────────
    figure(fig_offset + 4)
    clf
    subplot(2,1,1)
    plot(kT, Xs_k_mech(1,:), 'b-', 'LineWidth', 1.5), hold on
    plot(kT, Xs_k_mech(2,:), 'r-', 'LineWidth', 1.5)
    grid on
    ylabel('x_s(k)  (deg C)')
    title(['SS-MPC Mech – ' run_label ' – SS Target States x_s(k)'])
    legend('x_{s,1}(k)', 'x_{s,2}(k)', 'Location', 'Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(2,1,2)
    stairs(kT, Us_k_mech(1,:)', 'b-', 'LineWidth', 1.5), hold on
    stairs(kT, Us_k_mech(2,:)', 'r-', 'LineWidth', 1.5)
    grid on
    xlabel('Sample k')
    ylabel('u_s(k)  (% heater)')
    title(['SS-MPC Mech – ' run_label ' – SS Target Inputs u_s(k)'])
    legend('u_{s,1}(k)', 'u_{s,2}(k)', 'Location', 'Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── Save per-run results ──────────────────────────────────────────────
    if Noise_ON == 0
        save TCL_MPC_Mech_SimResults_NoiseFree ...
            kT t_sec Yk_mech Uk_mech Rk_f_mech ek_f_mech Xs_k_mech Us_k_mech ...
            ISE_T1 ISE_T2 IAE_T1 IAE_T2 Ts Ys Us Xs n_op n_ip N ...
            Np Nc Wx Wdelu beta_r alfa_e Plant_Sim
        fprintf('Saved: TCL_MPC_Mech_SimResults_NoiseFree.mat\n');
    else
        save TCL_MPC_Mech_SimResults_Noisy ...
            kT t_sec Yk_mech Uk_mech Rk_f_mech ek_f_mech Xs_k_mech Us_k_mech ...
            ISE_T1 ISE_T2 IAE_T1 IAE_T2 Ts Ys Us Xs n_op n_ip N ...
            Np Nc Wx Wdelu beta_r alfa_e Plant_Sim
        fprintf('Saved: TCL_MPC_Mech_SimResults_Noisy.mat\n');
    end

end  % Noise_ON loop

fprintf('\nAll SS-MPC simulations complete.\n');

% =========================================================================
%  END OF SCRIPT
% =========================================================================
