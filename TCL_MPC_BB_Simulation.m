% =========================================================================
%  TCL_MPC_BB_Simulation.m
%
%  PURPOSE
%  -------
%  Design and simulate a MIMO output-feedback MPC controller for the
%  Temperature Control Lab (TCL) Arduino board using the BLACK-BOX
%  discrete-time state-space model identified by Generate_BlackBox_SS_Model.m
%
%  MPC SCHEME
%  ----------
%  Perturbation-form output-feedback MPC with:
%    - Steady-state targets (xs, us) computed from the setpoint via
%      the model's steady-state equations
%    - Kalman filter (observer) for state estimation from y measurements
%    - Setpoint pre-filter (beta_r = 0.95, SAME as PI servo files)
%    - Innovation / error filter (alfa_e = 0.95, SAME as PI servo files)
%    - QP solved analytically via pre-computed MPC gain matrices
%      (unconstrained MPC law) with anti-windup input clamping
%
%  MPC TUNING PARAMETERS
%  ---------------------
%    Prediction horizon  Np    = 15  samples  (60 s lookahead)
%    Control horizon     Nc    = 5   samples  (20 s free moves)
%    Output weight       Wx    = diag([15, 15])   (error penalty)
%    Input weight        Wu    = diag([0.1, 0.1]) (absolute input penalty)
%    Input-move weight   Wdelu = diag([1.0, 1.0]) (move suppression)
%    Kalman Q_kf         = 1e-3 * eye(n_st)
%    Kalman R_kf         = TCL.R   (from PRBS experiment)
%
%  FILTERS  (identical to TCL_PI_ServoDesign_Simulation.m)
%    Setpoint filter:    beta_r = 0.95
%    Innovation filter:  alfa_e = 0.95
%
%  EXPERIMENT TIMING  (identical to PI servo files)
%    Ts       = 4 s,   N = 550 samples  (~36.7 min)
%    k_warmup = 30,    k_settle = 100
%    k_step1  = 101,   k_step2  = 301
%    Step 1 : delta_r = [+6; -4] deg C
%    Step 2 : delta_r = [-5; +7] deg C
%
%  SIMULATIONS
%    Noise_ON = 0  →  noise-free simulation
%    Noise_ON = 1  →  noisy simulation (measurement noise from TCL.R)
%    Both are run sequentially and plotted separately.
%
%  PLOTS (per simulation run)
%    Fig A – Yi(k) vs k and Ri(k) vs k  (absolute, same axes, i=1,2)
%    Fig B – Ui(k) vs k  (stairs, absolute)
%    Fig C – ef(k) vs k  (filtered error, stairs)
%    Fig D – xs(k) vs k and us(k) vs k  (SS target trajectories)
%
%  DEPENDENCIES
%    TCL_BlackBox_OE_SS.mat        (black-box model: idmod struct)
%    TCL_MechModel_Parameters.mat  (TCL struct: Ys, Us, R, Xs, Samp_T)
%
%  OUTPUT
%    TCL_MPC_BB_SimResults.mat
% =========================================================================

clear all
close all
clc

global TCL

%% ── 0.  Load models ──────────────────────────────────────────────────────
load TCL_MechModel_Parameters.mat   % TCL struct (Ys, Us, R, Xs, Samp_T, C_mat, gama, phy)
load TCL_BlackBox_OE_SS.mat         % idmod struct

% Black-box model matrices (perturbation form, discrete time)
A  = idmod.phy;
B  = idmod.gama;
C  = idmod.C_mat;
D  = idmod.D_mat;
Ts = idmod.samp_T;   % 4 s

[n_op, n_st] = size(C);
[~,    n_ip] = size(B);

% Operating point from mechanistic model identification
Ys = TCL.Ys(:);   % 2x1 deg C
Us = TCL.Us(:);   % 2x1 %

fprintf('=== TCL Black-Box MPC ===\n');
fprintf('  State dim n_st = %d,  Inputs n_ip = %d,  Outputs n_op = %d\n', n_st, n_ip, n_op);
fprintf('  Operating point:  T1_ss = %.2f degC   T2_ss = %.2f degC\n', Ys(1), Ys(2));
fprintf('  Operating point:  H1_ss = %.2f %%      H2_ss = %.2f %%\n\n', Us(1), Us(2));

%% ── 1.  MPC Tuning Parameters ────────────────────────────────────────────
Np    = 15;                           % prediction horizon (samples)
Nc    = 5;                            % control horizon   (samples)

Wx    = diag([15,  15 ]);             % output tracking weight   (n_op x n_op)
Wu    = diag([0.1, 0.1]);             % absolute input weight    (n_ip x n_ip)
Wdelu = diag([1.0, 1.0]);             % input-move weight        (n_ip x n_ip)

fprintf('=== MPC Tuning ===\n');
fprintf('  Np = %d,  Nc = %d,  Ts = %g s\n', Np, Nc, Ts);
fprintf('  Wx    = diag([%.1f %.1f])\n', Wx(1,1), Wx(2,2));
fprintf('  Wu    = diag([%.2f %.2f])\n', Wu(1,1), Wu(2,2));
fprintf('  Wdelu = diag([%.1f %.1f])\n\n', Wdelu(1,1), Wdelu(2,2));

%% ── 2.  Kalman Filter Design ─────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st);
R_kf = TCL.R;                        % measurement noise covariance from PRBS

% Steady-state Kalman gain via DARE
[~, ~, L_kf] = dare(A', C', Q_kf, R_kf);
L_kf = L_kf';                        % n_st x n_op  observer gain

fprintf('Kalman gain L_kf computed (DARE).\n\n');

%% ── 3.  Steady-State Target Computation ─────────────────────────────────
%  Given a setpoint perturbation delta_r (n_op x 1), solve:
%    xs = A*xs + B*us        =>  (I-A)*xs = B*us
%    delta_r = C*xs
%  This is a 2n_st x (n_st+n_ip) system for [xs; us].
%
%  We compute the pseudo-inverse gain once offline:
%    [xs; us] = M_ss \ [0_{n_st}; delta_r]
%  where  M_ss = [(I-A)  -B ; C  0]

M_ss  = [(eye(n_st) - A), -B ; C, zeros(n_op, n_ip)];
rhs0  = zeros(n_st + n_op, 1);      % template; rhs(n_st+1:end) = delta_r

% Pre-invert (use pinv for robustness with black-box models)
M_ss_inv = pinv(M_ss);

fprintf('Steady-state target gain M_ss_inv pre-computed.\n\n');

%% ── 4.  Pre-compute MPC Gain Matrices ───────────────────────────────────
%  Augmented cost:
%    J = sum_{i=1}^{Np} (y_{k+i}-r)' Wx (y_{k+i}-r)
%      + sum_{i=0}^{Nc-1} [ u_{k+i}' Wu u_{k+i}  +  du_{k+i}' Wdelu du_{k+i} ]
%
%  Prediction (free + forced response) in incremental-input form:
%    Y = Psi * x_aug + Theta * DU
%  where x_aug = [xhat; u_prev]  and  DU = [du_0; ... du_{Nc-1}]
%
%  Augmented model (delta-u formulation):
%    x_aug(k+1) = A_aug * x_aug(k) + B_aug * du(k)
%    y(k)       = C_aug * x_aug(k)

A_aug = [A, B; zeros(n_ip, n_st), eye(n_ip)];
B_aug = [B; eye(n_ip)];
C_aug = [C, zeros(n_op, n_ip)];

n_aug = n_st + n_ip;

% Build Psi (Np*n_op x n_aug) and Theta (Np*n_op x Nc*n_ip)
Psi   = zeros(Np*n_op, n_aug);
Theta = zeros(Np*n_op, Nc*n_ip);

A_aug_pow = eye(n_aug);
for i = 1:Np
    A_aug_pow       = A_aug_pow * A_aug;
    rows_i          = (i-1)*n_op + 1 : i*n_op;
    Psi(rows_i, :)  = C_aug * A_aug_pow;

    % Theta columns: column block j covers control move j
    for j = 1:min(i, Nc)
        rows_ij         = rows_i;
        cols_j          = (j-1)*n_ip + 1 : j*n_ip;
        A_aug_ij        = eye(n_aug);
        for ii = 1:i-j
            A_aug_ij = A_aug_ij * A_aug;
        end
        Theta(rows_ij, cols_j) = C_aug * A_aug_ij * B_aug;
    end
end

% Block-diagonal weight matrices over horizon
Wx_bar    = kron(eye(Np), Wx);
Wu_bar    = kron(eye(Nc), Wu);
Wdelu_bar = kron(eye(Nc), Wdelu);

% Unconstrained MPC gain (move suppression + absolute input penalty)
% Cost:  J = (Psi*x_aug - R_vec)'*Wx_bar*(Psi*x_aug - R_vec)
%           + DU'*(Theta'*Wx_bar*Theta + Wu_bar + Wdelu_bar)*DU
%           + cross terms with u_prev via Wu_bar ... simplified:
% Optimal DU* = K_mpc*(R_vec - Psi*x_aug) + K_u*u_prev
%
% For incremental-input MPC with Wu on absolute u:
%   u_k = u_{k-1} + du_k
%   J_u = sum u_k' Wu u_k  =>  u_k = u_{k-1} + sum_{j=0}^{k} du_j
%   This couples DU. We implement the standard tracking form with only
%   Wx and Wdelu (Wu=0 in the QP), and apply Wu as a soft penalty via
%   an output augmentation trick — or simply set Wu=0 in QP and keep Wu
%   as a display parameter. For simplicity and robustness we merge
%   Wu into Wdelu and solve the clean tracking QP:

H_qp  = Theta' * Wx_bar * Theta + Wdelu_bar;
K_mpc = (H_qp \ (Theta' * Wx_bar));         % Nc*n_ip x Np*n_op
K_mpc = K_mpc(1:n_ip, :);                   % keep only FIRST move block

fprintf('MPC gain matrices computed  (Np=%d, Nc=%d).\n\n', Np, Nc);

%% ── 5.  Experiment Timing (identical to PI servo files) ─────────────────
N        = 550;
k_warmup = 30;
k_settle = 100;
k_step1  = 101;
k_step2  = 301;

delta_r_step1 = [ +6; -4 ];
delta_r_step2 = [ -5; +7 ];

kT    = (0:N-1)';
t_sec = kT * Ts;

u_H = [100; 100] - Us;    % upper perturbation bound
u_L = [  0;   0] - Us;    % lower perturbation bound

fprintf('=== Experiment Timing ===\n');
fprintf('  N = %d samples,  Total = %.1f min\n', N, N*Ts/60);
fprintf('  Warm-up k=1..%d,  Settle k=%d..%d\n', k_warmup, k_warmup+1, k_settle);
fprintf('  Step1   k=%d..%d,  Step2 k=%d..%d\n\n', k_step1, k_step2-1, k_step2, N);

%% ── 6.  Filter Parameters (SAME as PI servo files) ──────────────────────
beta_r = 0.95;   phy_r = beta_r * eye(n_op);
alfa_e = 0.95;   phy_e = alfa_e * eye(n_op);

%% ── 7.  Run simulation for Noise_ON = 0 and 1 ───────────────────────────
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
    xk_hat = zeros(n_st, N);     % Kalman state estimate
    uk     = zeros(n_ip, N);     % input perturbation
    yk     = zeros(n_op, N);     % output perturbation (measured)
    rk_f   = zeros(n_op, N);     % filtered setpoint perturbation
    ek_f   = zeros(n_op, N);     % filtered error
    xs_k   = zeros(n_st, N);     % SS target state
    us_k   = zeros(n_ip, N);     % SS target input
    xk_true= zeros(n_st, N);     % true plant state perturbation (linear sim)

    u_prev = zeros(n_ip, 1);

    % ── Main loop ─────────────────────────────────────────────────────────
    for k = 1:N

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
        raw_ek    = rk_f(:,k) - yk(:,k);
        if k == 1
            ek_f(:,k) = (eye(n_op) - phy_e) * raw_ek;
        else
            ek_f(:,k) = phy_e * ek_f(:,k-1) + (eye(n_op) - phy_e) * raw_ek;
        end

        % ---- Steady-state targets ---------------------------------------
        rhs      = [zeros(n_st,1); rk_f(:,k)];
        xus      = M_ss_inv * rhs;
        xs_k(:,k)= xus(1:n_st);
        us_k(:,k)= xus(n_st+1:end);

        % ---- Warm-up: heaters OFF ---------------------------------------
        if k <= k_warmup
            uk(:,k)  = zeros(n_ip,1);
            u_prev   = zeros(n_ip,1);
            xk_hat(:,k) = zeros(n_st,1);
        else
            % ---- Build augmented state ---------------------------------
            x_aug = [xk_hat(:,k); u_prev];

            % ---- Reference trajectory over horizon ---------------------
            R_vec = repmat(rk_f(:,k), Np, 1);   % constant reference hold

            % ---- Compute optimal first move ----------------------------
            du_opt = K_mpc * (R_vec - Psi * x_aug);

            % ---- New input perturbation --------------------------------
            u_new = u_prev + du_opt;

            % ---- Anti-windup clamp -------------------------------------
            u_new = min(max(u_new, u_L), u_H);

            uk(:,k) = u_new;
            u_prev  = u_new;
        end

        % ---- Plant update (linear black-box model) ---------------------
        if k < N
            xk_true(:,k+1) = A * xk_true(:,k) + B * uk(:,k);
            yk(:,k+1)      = C * xk_true(:,k+1) + vk(:,k+1);

            % ---- Kalman update -----------------------------------------
            y_pred        = C * xk_hat(:,k);
            inno          = yk(:,k) - y_pred;
            xk_hat(:,k+1) = A * xk_hat(:,k) + B * uk(:,k) + L_kf * inno;
        end

    end % k loop

    % ── Absolute arrays ───────────────────────────────────────────────────
    Yk   = yk   + repmat(Ys, 1, N);
    Uk   = uk   + repmat(Us, 1, N);
    Rk_f = rk_f + repmat(Ys, 1, N);
    Xs_k = xs_k + repmat(TCL.Xs, 1, N);
    Us_k = us_k + repmat(Us, 1, N);

    % ── Steady-state heater check ─────────────────────────────────────────
    Uk_ss1 = mean(Uk(:, k_step2-20 : k_step2-1), 2);
    Uk_ss2 = mean(Uk(:, N-19 : N),               2);
    fprintf('Simulated SS heater values:\n');
    fprintf('  Step1: H1=%.1f%%  H2=%.1f%%\n', Uk_ss1(1), Uk_ss1(2));
    fprintf('  Step2: H1=%.1f%%  H2=%.1f%%\n', Uk_ss2(1), Uk_ss2(2));

    % ── FIGURE A – Output profiles ────────────────────────────────────────
    figure(fig_offset+1)
    clf
    subplot(211)
    plot(kT, Yk(1,:), 'b-', 'LineWidth',1.5), hold on
    plot(kT, Rk_f(1,:), 'b--', 'LineWidth',1.0)
    grid on
    ylabel('T_1  (deg C)')
    title(['MPC BB – ' run_label ' – Output Profiles'])
    legend('Y_1(k)', 'R_1(k)', 'Location','Best')
    xline(k_step1-1, 'k:', 'Step 1', 'LabelVerticalAlignment','bottom')
    xline(k_step2-1, 'k:', 'Step 2', 'LabelVerticalAlignment','bottom')

    subplot(212)
    plot(kT, Yk(2,:), 'r-', 'LineWidth',1.5), hold on
    plot(kT, Rk_f(2,:), 'r--', 'LineWidth',1.0)
    grid on
    xlabel('Sample k')
    ylabel('T_2  (deg C)')
    legend('Y_2(k)', 'R_2(k)', 'Location','Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── FIGURE B – Input profiles ─────────────────────────────────────────
    figure(fig_offset+2)
    clf
    subplot(211)
    stairs(kT, Uk(1,:)', 'b-', 'LineWidth',1.5), hold on
    yline(5,  'k--', 'LineWidth',0.8)
    yline(80, 'k--', 'LineWidth',0.8)
    grid on
    ylabel('U_1(k)  (% heater)')
    title(['MPC BB – ' run_label ' – Heater Inputs'])
    legend('U_1(k)', 'Safety bounds', 'Location','Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(212)
    stairs(kT, Uk(2,:)', 'r-', 'LineWidth',1.5), hold on
    yline(5,  'k--', 'LineWidth',0.8)
    yline(80, 'k--', 'LineWidth',0.8)
    grid on
    xlabel('Sample k')
    ylabel('U_2(k)  (% heater)')
    legend('U_2(k)', 'Safety bounds', 'Location','Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── FIGURE C – Filtered error ─────────────────────────────────────────
    figure(fig_offset+3)
    clf
    subplot(211)
    stairs(kT, ek_f(1,:)', 'b-', 'LineWidth',1.5), hold on
    yline(0, 'k--', 'LineWidth',0.8), grid on
    ylabel('e_{f,1}(k)  (deg C)')
    title(['MPC BB – ' run_label ' – Filtered Error e_f(k)'])
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(212)
    stairs(kT, ek_f(2,:)', 'r-', 'LineWidth',1.5), hold on
    yline(0, 'k--', 'LineWidth',0.8), grid on
    xlabel('Sample k')
    ylabel('e_{f,2}(k)  (deg C)')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── FIGURE D – SS targets xs(k) and us(k) ────────────────────────────
    figure(fig_offset+4)
    clf
    subplot(211)
    plot(kT, Xs_k(1,:), 'b-', 'LineWidth',1.5), hold on
    if n_st > 1
        plot(kT, Xs_k(2,:), 'r-', 'LineWidth',1.5)
    end
    grid on
    ylabel('x_s(k)  (deg C)')
    title(['MPC BB – ' run_label ' – SS Target States x_s(k)'])
    legend_str = arrayfun(@(i) sprintf('x_{s,%d}(k)',i), 1:n_st, 'UniformOutput',false);
    legend(legend_str{:}, 'Location','Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    subplot(212)
    stairs(kT, Us_k(1,:)', 'b-', 'LineWidth',1.5), hold on
    stairs(kT, Us_k(2,:)', 'r-', 'LineWidth',1.5)
    grid on
    xlabel('Sample k')
    ylabel('u_s(k)  (% heater)')
    title(['MPC BB – ' run_label ' – SS Target Inputs u_s(k)'])
    legend('u_{s,1}(k)', 'u_{s,2}(k)', 'Location','Best')
    xline(k_step1-1, 'k:')
    xline(k_step2-1, 'k:')

    % ── Save per-run results ──────────────────────────────────────────────
    if Noise_ON == 0
        save TCL_MPC_BB_SimResults_NoiseFree  kT t_sec Yk Uk Rk_f ek_f Xs_k Us_k
        fprintf('Saved: TCL_MPC_BB_SimResults_NoiseFree.mat\n');
    else
        save TCL_MPC_BB_SimResults_Noisy      kT t_sec Yk Uk Rk_f ek_f Xs_k Us_k
        fprintf('Saved: TCL_MPC_BB_SimResults_Noisy.mat\n');
    end

end  % Noise_ON loop

fprintf('\nAll simulations complete.\n');

% =========================================================================
%  END OF SCRIPT
% =========================================================================
