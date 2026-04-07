% =========================================================================
%  TCL_MPC_OT_Simulation.m
%
%  PURPOSE
%  -------
%  Design and simulate a MIMO OUTPUT-TRACKING MPC controller for the
%  Temperature Control Lab (TCL) using the MECHANISTIC linearised
%  discrete-time state-space model (TCL.phy, TCL.gama).
%
%  MPC SCHEME  –  Output-Tracking MPC (OT-MPC)
%  ---------------------------------------------
%  The QP cost penalises OUTPUT deviations from the setpoint reference
%  directly over the entire prediction horizon, plus input increments Du:
%
%    J = sum_{i=1}^{Np} (y(k+i) - r(k))' Wx (y(k+i) - r(k))
%      + sum_{j=0}^{Nc-1} Du(k+j)' Wdelu Du(k+j)
%
%  Prediction model (perturbation form, delta-u):
%    x_aug(k+1) = A_aug * x_aug(k) + B_aug * Du(k)
%    y(k+i)     = C_aug * x_aug(k+i)
%  where  x_aug = [x_hat ; u_prev]
%
%  Reference vector over horizon:
%    R_vec = [rk_f; rk_f; ... (Np times)]   in OUTPUT space (Np*n_op x 1)
%  The controller directly minimises predicted OUTPUT error vs setpoint.
%
%  Key difference vs State-Space MPC:
%    OT-MPC  : reference = rk_f  (output setpoint, repeated over horizon)
%    SS-MPC  : reference = xs_k  (state target, from SS equations)
%
%  CONTROLLER STRUCTURE
%  --------------------
%    - Augmented state  x_aug(k) = [x_hat(k); u(k-1)]
%    - Psi   : maps x_aug to predicted OUTPUT trajectory  (Np*n_op x n_aug)
%    - Theta : maps Delta_U to predicted OUTPUT trajectory (Np*n_op x Nc*n_ip)
%    - Unconstrained MPC gain K_mpc_ot from output cost
%    - Kalman filter for state estimation (DARE)
%    - SS targets (xs_k, us_k) computed for diagnostics / plotting only
%    - Anti-windup input clamping to [0, 100]%
%
%  MPC TUNING PARAMETERS
%  ---------------------
%    Prediction horizon  Np    = 15  samples  (60 s lookahead)
%    Control horizon     Nc    = 5   samples  (20 s free moves)
%    Output weight       Wx    = diag([15, 15])   (output error penalty)
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
%    Fig 1/5 – Yi(k) vs k  and  Ri(k) vs k
%    Fig 2/6 – Ui(k) vs k  (heater inputs, stairs)
%    Fig 3/7 – ef(k) vs k  (filtered error)
%    Fig 4/8 – xs(k) vs k  and  us(k) vs k  (SS target trajectories)
%
%  DEPENDENCIES
%    TCL_MechModel_Parameters.mat  (TCL struct)
%    TCL_Dynamics.m                (nonlinear ODE, needed if Plant_Sim=0)
%
%  OUTPUT
%    TCL_MPC_OT_SimResults_NoiseFree.mat
%    TCL_MPC_OT_SimResults_Noisy.mat
%      Variables (suffix _ot): Yk_ot, Uk_ot, Rk_f_ot,
%                              ek_f_ot, Xs_k_ot, Us_k_ot
% =========================================================================

clear all
close all
clc

global TCL

%% ── 0.  Load mechanistic model parameters ────────────────────────────────
load TCL_MechModel_Parameters.mat    % TCL struct

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

fprintf('=== TCL OUTPUT-TRACKING MPC (Mechanistic Model) ===\n');
fprintf('  MPC scheme : OT-MPC  (output error cost, rk_f reference over horizon)\n');
fprintf('  n_st=%d  n_ip=%d  n_op=%d   Ts=%g s\n', n_st, n_ip, n_op, Ts);
fprintf('  Op. point:  T1_ss=%.2f degC   T2_ss=%.2f degC\n', Ys(1), Ys(2));
fprintf('  Op. point:  H1_ss=%.2f %%     H2_ss=%.2f %%\n\n', Us(1), Us(2));

%% ── 1.  MPC Tuning Parameters ────────────────────────────────────────────
Np    = 15;                       % prediction horizon (samples)
Nc    = 5;                        % control horizon   (samples)
Wx    = diag([15,  15 ]);         % OUTPUT error weight  (n_op x n_op)
Wdelu = diag([1.0, 1.0]);         % input-move weight    (n_ip x n_ip)

fprintf('=== MPC Tuning ===\n');
fprintf('  Np=%d  Nc=%d  Ts=%g s\n', Np, Nc, Ts);
fprintf('  Wx (output weight)   = diag([%.1f %.1f])\n', Wx(1,1), Wx(2,2));
fprintf('  Wdelu (move weight)  = diag([%.1f %.1f])\n\n', Wdelu(1,1), Wdelu(2,2));

%% ── 2.  Kalman Filter Design ─────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st);
R_kf = TCL.R;                    % measurement noise covariance (from PRBS)

[~, ~, L_kf] = dare(A', C', Q_kf, R_kf);
L_kf = L_kf';                   % n_st x n_op  observer gain

fprintf('Kalman gain L_kf computed (DARE).\n\n');

%% ── 3.  Steady-State Target Computation (diagnostics / plotting only) ────
%  In OT-MPC, xs_k is NOT used in the QP reference.
%  The controller tracks rk_f directly through the OUTPUT prediction.
%  xs_k and us_k are retained only for comparison plots.
M_ss     = [(eye(n_st) - A), -B ; C, zeros(n_op, n_ip)];
M_ss_inv = pinv(M_ss);

fprintf('SS target solver M_ss_inv computed (used for diagnostics only).\n\n');

%% ── 4.  Augmented State-Space (delta-u formulation) ─────────────────────
%  x_aug(k+1) = A_aug * x_aug(k) + B_aug * Du(k)
%  y(k)       = C_aug * x_aug(k)
%  where x_aug = [x_hat(k); u(k-1)],  Du = u(k) - u(k-1)

n_aug = n_st + n_ip;
A_aug = [A,                 B             ;
         zeros(n_ip, n_st), eye(n_ip)     ];
B_aug = [B        ;
         eye(n_ip)];

% C_aug maps x_aug -> OUTPUTS (n_op x n_aug)
C_aug = [C,  zeros(n_op, n_ip)];

%% ── 5.  Build Prediction Matrices (OUTPUT-based) ─────────────────────────
%  Predicted output trajectory:
%    Y_pred = Psi * x_aug(k) + Theta * Delta_U
%  where Y_pred = [y(k+1); y(k+2); ...; y(k+Np)]   size Np*n_op x 1
%
%  This is the OUTPUT prediction — the reference R_vec is also in output
%  space: R_vec = [rk_f; rk_f; ...(Np)] size Np*n_op x 1

Psi   = zeros(Np * n_op,  n_aug);
Theta = zeros(Np * n_op,  Nc * n_ip);

A_aug_pow = eye(n_aug);
for i = 1 : Np
    A_aug_pow     = A_aug_pow * A_aug;
    rows_i        = (i-1)*n_op + 1 : i*n_op;
    Psi(rows_i,:) = C_aug * A_aug_pow;       % output prediction rows

    for j = 1 : min(i, Nc)
        cols_j = (j-1)*n_ip + 1 : j*n_ip;
        A_ij   = eye(n_aug);
        for ii = 1 : i-j
            A_ij = A_ij * A_aug;
        end
        Theta(rows_i, cols_j) = C_aug * A_ij * B_aug;
    end
end

%% ── 6.  Pre-compute OT-MPC Gain K_mpc_ot ────────────────────────────────
%  Output cost over horizon:
%    J_y  = (Y_pred - R_vec)' Wx_bar (Y_pred - R_vec)
%         = DU' Theta' Wx_bar Theta DU
%           + 2*(Psi*x_aug - R_vec)' Wx_bar Theta DU + const
%
%  Unconstrained optimum:
%    DU* = (Theta' Wx_bar Theta + Wdelu_bar)^{-1}
%           * Theta' Wx_bar * (R_vec - Psi * x_aug)
%
%  Only the FIRST block row DU*(1) is applied (receding horizon).

Wx_bar    = kron(eye(Np), Wx);        % Np*n_op x Np*n_op
Wdelu_bar = kron(eye(Nc), Wdelu);     % Nc*n_ip x Nc*n_ip

H_ot      = Theta' * Wx_bar * Theta + Wdelu_bar;
H_ot      = (H_ot + H_ot') / 2;              % enforce symmetry
K_mpc_ot  = H_ot \ (Theta' * Wx_bar);        % Nc*n_ip x Np*n_op
K_mpc_ot  = K_mpc_ot(1:n_ip, :);             % first move block [n_ip x Np*n_op]

fprintf('OT-MPC gain K_mpc_ot computed  (size %dx%d).\n', ...
        size(K_mpc_ot,1), size(K_mpc_ot,2));
fprintf('  Reference vector : rk_f (output setpoint) repeated Np times.\n\n');

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
for Noise_ON = 0 : 1

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

    % ── Allocate storage ─────────────────────────────────────────────────
    xk_hat  = zeros(n_st, N);
    xk_true = zeros(n_st, N);
    uk      = zeros(n_ip, N);
    yk      = zeros(n_op, N);
    rk_f    = zeros(n_op, N);
    ek_f    = zeros(n_op, N);
    xs_k    = zeros(n_st, N);    % for diagnostics only
    us_k    = zeros(n_ip, N);    % for diagnostics only

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

        % ---- SS targets (diagnostics / plotting only) -------------------
        rhs       = [zeros(n_st, 1); rk_f(:,k)];
        xus       = M_ss_inv * rhs;
        xs_k(:,k) = xus(1 : n_st);
        us_k(:,k) = xus(n_st+1 : end);

        % ---- Warm-up: heaters OFF ----------------------------------------
        if k <= k_warmup
            uk(:,k)     = zeros(n_ip, 1);
            u_prev      = zeros(n_ip, 1);
            xk_hat(:,k) = zeros(n_st, 1);

        else
            % ---- Build augmented state ----------------------------------
            x_aug = [xk_hat(:,k); u_prev];

            % ---- OT-MPC: reference = rk_f repeated Np times in OUTPUT space
            %  R_vec is [rk_f; rk_f; ...(Np)] in OUTPUT space (Np*n_op x 1)
            %  This is the defining difference from SS-MPC (which uses xs_k)
            R_vec = repmat(rk_f(:,k), Np, 1);

            % ---- Optimal first move (OT-MPC law) ------------------------
            %  Du* = K_mpc_ot * (R_vec - Psi * x_aug)
            du_opt  = K_mpc_ot * (R_vec - Psi * x_aug);
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
    Yk_ot   = yk   + repmat(Ys, 1, N);
    Uk_ot   = uk   + repmat(Us, 1, N);
    Rk_f_ot = rk_f + repmat(Ys, 1, N);
    ek_f_ot = ek_f;
    Xs_k_ot = xs_k + repmat(Xs, 1, N);
    Us_k_ot = us_k + repmat(Us, 1, N);

    % ── Performance metrics ───────────────────────────────────────────────
    ISE_T1 = sum((Yk_ot(1,:) - Rk_f_ot(1,:)).^2) * Ts;
    ISE_T2 = sum((Yk_ot(2,:) - Rk_f_ot(2,:)).^2) * Ts;
    IAE_T1 = sum(abs(Yk_ot(1,:) - Rk_f_ot(1,:))) * Ts;
    IAE_T2 = sum(abs(Yk_ot(2,:) - Rk_f_ot(2,:))) * Ts;

    fprintf('=== Tracking Performance [%s] ===\n', run_label);
    fprintf('  ISE T1 = %8.2f   |   ISE T2 = %8.2f\n', ISE_T1, ISE_T2);
    fprintf('  IAE T1 = %8.2f   |   IAE T2 = %8.2f\n', IAE_T1, IAE_T2);

    % ── FIGURE 1 – Output profiles ────────────────────────────────────────
    figure(fig_offset + 1)
    clf
    subplot(2,1,1)
    plot(kT, Yk_ot(1,:), 'b-', 'LineWidth', 1.5), hold on
    plot(kT, Rk_f_ot(1,:), 'b--', 'LineWidth', 1.0)
    grid on
    ylabel('T_1  (deg C)')
    title(['OT-MPC Mech – ' run_label ' – Output Profiles'])
    legend('Y_1(k)', 'R_1(k)', 'Location', 'Best')
    xline(k_step1-1, 'k:', 'Step 1', 'LabelVerticalAlignment', 'bottom')
    xline(k_step2-1, 'k:', 'Step 2', 'LabelVerticalAlignment', 'bottom')

    subplot(2,1,2)
    plot(kT, Yk_ot(2,:), 'r-', 'LineWidth', 1.5), hold on
    plot(kT, Rk_f_ot(2,:), 'r--', 'LineWidth', 1.0)
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
    stairs(kT, Uk_ot(1,:)', 'b-', 'LineWidth', 1.5), hold on
    yline(5,  'k--', 'LineWidth', 0.8)
    yline(80, 'k--', 'LineWidth', 0.8)
    grid on
    ylabel('U_1(k)  (% heater)')
    title(['OT-MPC Mech – ' run_label ' – Heater Inputs'])
    legend('U_1(k)', 'Safety bounds', 'Location', 'Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(2,1,2)
    stairs(kT, Uk_ot(2,:)', 'r-', 'LineWidth', 1.5), hold on
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
    stairs(kT, ek_f_ot(1,:)', 'b-', 'LineWidth', 1.5), hold on
    yline(0, 'k--', 'LineWidth', 0.8), grid on
    ylabel('e_{f,1}(k)  (deg C)')
    title(['OT-MPC Mech – ' run_label ' – Filtered Error e_f(k)'])
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(2,1,2)
    stairs(kT, ek_f_ot(2,:)', 'r-', 'LineWidth', 1.5), hold on
    yline(0, 'k--', 'LineWidth', 0.8), grid on
    xlabel('Sample k')
    ylabel('e_{f,2}(k)  (deg C)')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── FIGURE 4 – SS targets xs(k) and us(k) ────────────────────────────
    figure(fig_offset + 4)
    clf
    subplot(2,1,1)
    plot(kT, Xs_k_ot(1,:), 'b-', 'LineWidth', 1.5), hold on
    plot(kT, Xs_k_ot(2,:), 'r-', 'LineWidth', 1.5)
    grid on
    ylabel('x_s(k)  (deg C)')
    title(['OT-MPC Mech – ' run_label ' – SS Target States x_s(k)  [diagnostics]'])
    legend('x_{s,1}(k)', 'x_{s,2}(k)', 'Location', 'Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(2,1,2)
    stairs(kT, Us_k_ot(1,:)', 'b-', 'LineWidth', 1.5), hold on
    stairs(kT, Us_k_ot(2,:)', 'r-', 'LineWidth', 1.5)
    grid on
    xlabel('Sample k')
    ylabel('u_s(k)  (% heater)')
    title(['OT-MPC Mech – ' run_label ' – SS Target Inputs u_s(k)  [diagnostics]'])
    legend('u_{s,1}(k)', 'u_{s,2}(k)', 'Location', 'Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── Save per-run results ──────────────────────────────────────────────
    if Noise_ON == 0
        save TCL_MPC_OT_SimResults_NoiseFree ...
            kT t_sec Yk_ot Uk_ot Rk_f_ot ek_f_ot Xs_k_ot Us_k_ot ...
            ISE_T1 ISE_T2 IAE_T1 IAE_T2 Ts Ys Us Xs n_op n_ip N ...
            Np Nc Wx Wdelu beta_r alfa_e Plant_Sim
        fprintf('Saved: TCL_MPC_OT_SimResults_NoiseFree.mat\n');
    else
        save TCL_MPC_OT_SimResults_Noisy ...
            kT t_sec Yk_ot Uk_ot Rk_f_ot ek_f_ot Xs_k_ot Us_k_ot ...
            ISE_T1 ISE_T2 IAE_T1 IAE_T2 Ts Ys Us Xs n_op n_ip N ...
            Np Nc Wx Wdelu beta_r alfa_e Plant_Sim
        fprintf('Saved: TCL_MPC_OT_SimResults_Noisy.mat\n');
    end

end  % Noise_ON loop

fprintf('\nAll OT-MPC simulations complete.\n');

% =========================================================================
%  END OF SCRIPT
% =========================================================================
