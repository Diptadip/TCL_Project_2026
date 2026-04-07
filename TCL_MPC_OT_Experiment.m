% =========================================================================
%  TCL_MPC_OT_Experiment.m
%
%  PURPOSE
%  -------
%  Runs the Output-Tracking MPC servo-control experiment on the TCL
%  Arduino board using the MECHANISTIC discrete-time state-space model
%  (TCL.phy, TCL.gama) obtained from linearization.
%
%  MPC SCHEME
%  ----------
%  Output-feedback MPC (augmented-state / delta-u formulation):
%    - Kalman filter reconstructs states from T1, T2 measurements
%    - Steady-state targets (xs, us) computed online from filtered setpoint
%    - Unconstrained MPC law (pre-computed gain K_mpc) + anti-windup clamp
%    - Setpoint pre-filter  beta_r = 0.95  (SAME as PI & BB-MPC experiments)
%    - Innovation filter    alfa_e = 0.95  (SAME as PI & BB-MPC experiments)
%
%  SETPOINT SCHEDULE  (identical to all project servo files)
%    Phase 1  k = 1   … 30   Warm-up: heaters OFF
%    Phase 2  k = 31  … 100  MPC active, setpoint = Ys  (settle)
%    Phase 3  k = 101 … 300  Step 1: T1 +6 deg C,  T2 -4 deg C
%    Phase 4  k = 301 … 550  Step 2: T1 -5 deg C,  T2 +7 deg C
%
%  Total duration: 550 × 4 s = 2200 s ≈ 36.7 min
%
%  MPC PARAMETERS  (same as TCL_MPC_OT_Simulation.m)
%    Np = 15,  Nc = 5
%    Wx = diag([15, 15]),  Wdelu = diag([1.0, 1.0])
%    Kalman: Q_kf = 1e-3*I,  R_kf = TCL.R
%
%  PLOTS
%    Fig 1 – Live strip chart (updated every 15 samples = 60 s)
%    Fig 2 – Yi(k) vs k  and  Ri(k) vs k   (post, i=1,2)
%    Fig 3 – Ui(k) vs k  (stairs, post)
%    Fig 4 – ef(k) vs k  (stairs, post)
%    Fig 5 – xs(k) vs k  and  us(k) vs k   (post)
%
%  DEPENDENCIES
%    tclab_N.m
%    TCL_MechModel_Parameters.mat  (must contain TCL.phy, TCL.gama, TCL.R)
%    TCL_MPC_OT_Simulation.m       (run first to validate gains offline)
%
%  OUTPUT
%    TCL_MPC_OT_ExpResults.mat
%      Variables (suffix _ot) compatible with TCL_Compare_Experiment.m:
%        Yk_ot, Uk_ot, Rk_f_ot, ek_f_ot, Xs_k_ot, Us_k_ot, kT, Ts
% =========================================================================

close all; clear all; clc

%% ── 0.  Connect to Arduino ───────────────────────────────────────────────
tclab_N;

%% ── 1.  Load mechanistic model ───────────────────────────────────────────
disp('Loading mechanistic model and computing MPC parameters ...')
load TCL_MechModel_Parameters.mat    % TCL struct

A  = TCL.phy;                        % [n_st x n_st]
B  = TCL.gama;                       % [n_st x n_ip]
C  = TCL.C_mat;                      % [n_op x n_st]
D  = TCL.D_mat;                      % [n_op x n_ip]
Ts = TCL.Samp_T;                     % 4 s

[n_op, n_st] = size(C);
[~,    n_ip] = size(B);

Ys = TCL.Ys(:);                      % steady-state temperatures [deg C]
Us = TCL.Us(:);                      % steady-state heater inputs [%]
Xs = TCL.Xs(:);                      % steady-state states       [deg C]

fprintf('Operating point:  T1_ss=%.2f degC   T2_ss=%.2f degC\n', Ys(1), Ys(2));
fprintf('                  H1_ss=%.2f %%       H2_ss=%.2f %%\n\n', Us(1), Us(2));

%% ── 2.  MPC Tuning Parameters ────────────────────────────────────────────
Np    = 15;
Nc    = 5;
Wx    = diag([15,  15 ]);
Wdelu = diag([1.0, 1.0]);

%% ── 3.  Kalman Filter ────────────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st);
R_kf = TCL.R;
[~, ~, L_kf] = dare(A', C', Q_kf, R_kf);
L_kf = L_kf';                        % [n_st x n_op]

%% ── 4.  Steady-State Target Gain ─────────────────────────────────────────
M_ss     = [(eye(n_st) - A),  -B ;
             C,                zeros(n_op, n_ip)];
M_ss_inv = pinv(M_ss);

%% ── 5.  Pre-compute MPC Gain K_mpc ──────────────────────────────────────
%  Augmented state  xi(k) = [ x_hat(k) ; u(k-1) ]
n_aug = n_st + n_ip;
A_aug = [A,                 B              ;
         zeros(n_ip, n_st), eye(n_ip)      ];
B_aug = [B        ;
         eye(n_ip)];
C_aug = [C,  zeros(n_op, n_ip)];

Psi   = zeros(Np * n_op,  n_aug);
Theta = zeros(Np * n_op,  Nc * n_ip);

A_aug_pow = eye(n_aug);
for i = 1 : Np
    A_aug_pow     = A_aug_pow * A_aug;
    rows_i        = (i-1)*n_op + 1 : i*n_op;
    Psi(rows_i,:) = C_aug * A_aug_pow;
    for j = 1 : min(i, Nc)
        cols_j = (j-1)*n_ip + 1 : j*n_ip;
        A_ij   = eye(n_aug);
        for ii = 1 : i-j
            A_ij = A_ij * A_aug;
        end
        Theta(rows_i, cols_j) = C_aug * A_ij * B_aug;
    end
end

Wx_bar    = kron(eye(Np), Wx);
Wdelu_bar = kron(eye(Nc), Wdelu);
H_qp      = Theta' * Wx_bar * Theta + Wdelu_bar;
H_qp      = (H_qp + H_qp') / 2;     % enforce symmetry
K_mpc_full = H_qp \ (Theta' * Wx_bar);
K_mpc      = K_mpc_full(1:n_ip, :); % first move block only  [n_ip x Np*n_op]

disp('MPC gain K_mpc pre-computed.  Ready to start experiment.')

%% ── 6.  Experiment Timing ────────────────────────────────────────────────
N_samp   = 550;
k_warmup = 30;
k_settle = 100;
k_step1  = 101;
k_step2  = 301;

delta_r_step1 = [ +6; -4];
delta_r_step2 = [ -5; +7];

u_H = [100; 100] - Us;               % upper input bound (perturbation)
u_L = [  0;   0] - Us;               % lower input bound (perturbation)

fprintf('=== Experiment Setpoints ===\n');
fprintf('  Step 1: T1->%.1f degC  T2->%.1f degC\n', ...
        Ys(1)+delta_r_step1(1), Ys(2)+delta_r_step1(2));
fprintf('  Step 2: T1->%.1f degC  T2->%.1f degC\n\n', ...
        Ys(1)+delta_r_step2(1), Ys(2)+delta_r_step2(2));

%% ── 7.  Filter Parameters (same as all project servo files) ──────────────
beta_r   = 0.95;   phy_r = beta_r * eye(n_op);
alfa_e   = 0.95;   phy_e = alfa_e * eye(n_op);
rk_f_now = zeros(n_op, 1);          % running filtered setpoint perturbation
ek_f_now = zeros(n_op, 1);          % running filtered error

%% ── 8.  MPC and Kalman state memory ─────────────────────────────────────
xk_hat = zeros(n_st, 1);            % Kalman state estimate (perturbation)
u_prev = zeros(n_ip, 1);            % previous input perturbation

%% ── 9.  Data logging ────────────────────────────────────────────────────
t1s  = [];   t2s  = [];
h1s  = [];   h2s  = [];
R1s  = [];   R2s  = [];
e1s  = [];   e2s  = [];
xs1s = [];   xs2s = [];
us1s = [];   us2s = [];

%% ── 10.  Live figure ─────────────────────────────────────────────────────
figure(1)

%% ── 11.  Main real-time loop ─────────────────────────────────────────────
fprintf('Starting MPC-OT experiment.  Total: %.1f min\n', N_samp*Ts/60);
fprintf('Phase 1 (warm-up): heaters OFF for %d s\n\n', k_warmup*Ts);

for k = 1 : N_samp
    tic

    % ---- Read hardware temperatures -------------------------------------
    Tk    = zeros(n_op, 1);
    Tk(1) = T1C();
    Tk(2) = T2C();
    yk_meas = Tk - Ys;               % output perturbation from steady state

    % ---- Setpoint schedule ----------------------------------------------
    if k < k_step1
        setpt = zeros(n_op, 1);
    elseif k < k_step2
        setpt = delta_r_step1;
    else
        setpt = delta_r_step2;
    end

    % ---- Setpoint filter ------------------------------------------------
    rk_f_now = phy_r * rk_f_now + (eye(n_op) - phy_r) * setpt;

    % ---- Innovation / error filter (logged as ef) -----------------------
    raw_ek   = rk_f_now - yk_meas;
    ek_f_now = phy_e * ek_f_now + (eye(n_op) - phy_e) * raw_ek;

    % ---- Kalman correction (measurement update) -------------------------
    y_pred = C * xk_hat;
    inno   = yk_meas - y_pred;
    xk_hat = xk_hat + L_kf * inno;  % corrected state estimate

    % ---- Steady-state targets -------------------------------------------
    rhs  = [zeros(n_st, 1); rk_f_now];
    xus  = M_ss_inv * rhs;
    xs_t = xus(1 : n_st);
    us_t = xus(n_st+1 : end);

    % ---- Phase 1: warm-up with heaters OFF ------------------------------
    if k <= k_warmup
        uk_pert  = zeros(n_ip, 1);
        u_prev   = zeros(n_ip, 1);
        xk_hat   = zeros(n_st,  1);
        ek_f_now = zeros(n_op, 1);
        rk_f_now = zeros(n_op, 1);

    else
        % ---- MPC law (unconstrained, pre-computed gain) -----------------
        x_aug  = [xk_hat; u_prev];
        R_vec  = repmat(rk_f_now, Np, 1);   % constant-reference hold
        du_opt = K_mpc * (R_vec - Psi * x_aug);

        % ---- Anti-windup clamp  (0% <= H1,H2 <= 100%) ------------------
        u_new  = u_prev + du_opt;
        u_new  = min(max(u_new, u_L), u_H);

        uk_pert = u_new;
        u_prev  = u_new;
    end

    % ---- Absolute heater values to send to hardware ---------------------
    Uk_abs = uk_pert + Us;
    Rk_abs = rk_f_now + Ys;

    % ---- Send to hardware -----------------------------------------------
    h1(Uk_abs(1));
    h2(Uk_abs(2));

    % ---- Kalman time update (prediction for next step) ------------------
    xk_hat = A * xk_hat + B * uk_pert;

    % ---- Log sample data ------------------------------------------------
    t1s  = [t1s,  Tk(1)];
    t2s  = [t2s,  Tk(2)];
    h1s  = [h1s,  Uk_abs(1)];
    h2s  = [h2s,  Uk_abs(2)];
    R1s  = [R1s,  Rk_abs(1)];
    R2s  = [R2s,  Rk_abs(2)];
    e1s  = [e1s,  ek_f_now(1)];
    e2s  = [e2s,  ek_f_now(2)];
    xs1s = [xs1s, xs_t(1) + Xs(1)];   % absolute SS target state
    xs2s = [xs2s, xs_t(2) + Xs(2)];
    us1s = [us1s, us_t(1) + Us(1)];   % absolute SS target input
    us2s = [us2s, us_t(2) + Us(2)];

    % ---- Phase banners in console ---------------------------------------
    if k == k_warmup + 1
        fprintf('Phase 2: MPC active, settling at Ys ...\n');
    elseif k == k_step1
        fprintf('Phase 3: Step 1  (T1%+.0f degC, T2%+.0f degC)\n', ...
                delta_r_step1(1), delta_r_step1(2));
    elseif k == k_step2
        fprintf('Phase 4: Step 2  (T1%+.0f degC, T2%+.0f degC)\n', ...
                delta_r_step2(1), delta_r_step2(2));
    end

    % ---- Live strip chart (every 15 samples = 60 s) ---------------------
    if rem(k, 15) == 0
        n_now = length(t1s);
        kk    = 0 : n_now-1;

        clf(figure(1))
        subplot(2,1,1)
        plot(kk, t1s, 'r-',  kk, R1s, 'r--', ...
             kk, t2s, 'b-',  kk, R2s, 'b--'), grid on
        ylabel('Temperature (deg C)')
        legend('Y_1','R_1','Y_2','R_2','Location','NorthWest')
        title(sprintf('TCL MPC-OT Experiment  –  k = %d / %d', k, N_samp))

        subplot(2,1,2)
        stairs(kk, h1s, 'r-', 'LineWidth', 2), hold on
        stairs(kk, h2s, 'b-', 'LineWidth', 2), grid on
        yline(  0, 'k--', 'LineWidth', 0.8)
        yline(100, 'k--', 'LineWidth', 0.8)
        ylabel('Heater (%)')
        xlabel('Sample k')
        legend('H_1','H_2','Location','NorthWest')
        drawnow
    end

    % ---- Console output -------------------------------------------------
    elapsed = toc;
    fprintf('k=%3d | T1=%.2f T2=%.2f | H1=%.1f%% H2=%.1f%% | slack=%.3fs\n', ...
            k, Tk(1), Tk(2), Uk_abs(1), Uk_abs(2), Ts - elapsed);
    pause(max(0.01, Ts - elapsed))

end  % main real-time loop

%% ── 12.  Turn heaters OFF safely ─────────────────────────────────────────
h1(0);   h2(0);
fprintf('\nExperiment complete. Heaters turned OFF.\n\n');

%% ── 13.  Assemble Post-Experiment Arrays ─────────────────────────────────
N   = N_samp;
kT  = (0 : N-1)';

% Absolute temperature and reference  [1 x N]
Yk_ot(1,:)   = t1s;
Yk_ot(2,:)   = t2s;
Rk_f_ot(1,:) = R1s;
Rk_f_ot(2,:) = R2s;

% Absolute heater inputs  [1 x N]
Uk_ot(1,:)   = h1s;
Uk_ot(2,:)   = h2s;

% Filtered error  [1 x N]  (perturbation)
ek_f_ot(1,:) = e1s;
ek_f_ot(2,:) = e2s;

% Steady-state target trajectories  (absolute)
Xs_k_ot(1,:) = xs1s;
Xs_k_ot(2,:) = xs2s;
Us_k_ot(1,:) = us1s;
Us_k_ot(2,:) = us2s;

%% ── 14.  Performance Metrics ─────────────────────────────────────────────
ISE_T1 = sum((Yk_ot(1,:) - Rk_f_ot(1,:)).^2) * Ts;
ISE_T2 = sum((Yk_ot(2,:) - Rk_f_ot(2,:)).^2) * Ts;
IAE_T1 = sum(abs(Yk_ot(1,:) - Rk_f_ot(1,:))) * Ts;
IAE_T2 = sum(abs(Yk_ot(2,:) - Rk_f_ot(2,:))) * Ts;

fprintf('=== Tracking Performance ===\n');
fprintf('  ISE T1 = %8.2f   |   ISE T2 = %8.2f\n', ISE_T1, ISE_T2);
fprintf('  IAE T1 = %8.2f   |   IAE T2 = %8.2f\n\n', IAE_T1, IAE_T2);

%% ── 15.  Post-Experiment Plots ────────────────────────────────────────────

% ── Figure 2 : Yi(k) vs k  and  Ri(k) vs k  (i = 1, 2) ──────────────────
figure('Name', 'MPC-OT Exp: Output Tracking', 'NumberTitle', 'off');

subplot(2,1,1)
plot(kT, Yk_ot(1,:), 'b-', 'LineWidth', 1.8), hold on
stairs(kT, Rk_f_ot(1,:), 'r--', 'LineWidth', 1.5)
grid on
ylabel('T_1 (°C)')
title(sprintf('MPC-OT Experiment — T_1   (ISE=%.1f, IAE=%.1f)', ISE_T1, IAE_T1))
legend('Y_1(k)','R_1(k)','Location','best')
xlim([0, N-1])

subplot(2,1,2)
plot(kT, Yk_ot(2,:), 'b-', 'LineWidth', 1.8), hold on
stairs(kT, Rk_f_ot(2,:), 'r--', 'LineWidth', 1.5)
grid on
ylabel('T_2 (°C)')
xlabel('Sample k')
title(sprintf('MPC-OT Experiment — T_2   (ISE=%.1f, IAE=%.1f)', ISE_T2, IAE_T2))
legend('Y_2(k)','R_2(k)','Location','best')
xlim([0, N-1])

% ── Figure 3 : Ui(k) vs k  (stairs, absolute heater %) ───────────────────
figure('Name', 'MPC-OT Exp: Heater Inputs', 'NumberTitle', 'off');

subplot(2,1,1)
stairs(kT, Uk_ot(1,:), 'k-', 'LineWidth', 1.5), hold on, grid on
yline(100, 'r--', 'H_{max}', 'LineWidth', 1.2)
yline(  0, 'b--', 'H_{min}', 'LineWidth', 1.2)
ylabel('U_1(k)  (%)'), title('Heater 1 Input — U_1(k)')
ylim([-5, 105]), xlim([0, N-1])

subplot(2,1,2)
stairs(kT, Uk_ot(2,:), 'k-', 'LineWidth', 1.5), hold on, grid on
yline(100, 'r--', 'H_{max}', 'LineWidth', 1.2)
yline(  0, 'b--', 'H_{min}', 'LineWidth', 1.2)
ylabel('U_2(k)  (%)'), xlabel('Sample k'), title('Heater 2 Input — U_2(k)')
ylim([-5, 105]), xlim([0, N-1])

% ── Figure 4 : ef(k) vs k  (filtered error, stairs) ──────────────────────
figure('Name', 'MPC-OT Exp: Filtered Error', 'NumberTitle', 'off');

subplot(2,1,1)
stairs(kT, ek_f_ot(1,:), 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5), hold on
yline(0, 'k--', 'LineWidth', 1.0), grid on
ylabel('e_{f,1}(k)  (°C)'), title('Filtered Tracking Error — e_{f,1}(k)')
xlim([0, N-1])

subplot(2,1,2)
stairs(kT, ek_f_ot(2,:), 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5), hold on
yline(0, 'k--', 'LineWidth', 1.0), grid on
ylabel('e_{f,2}(k)  (°C)'), xlabel('Sample k')
title('Filtered Tracking Error — e_{f,2}(k)')
xlim([0, N-1])

% ── Figure 5 : xs(k) vs k  and  us(k) vs k  (SS target trajectories) ─────
figure('Name', 'MPC-OT Exp: Steady-State Targets', 'NumberTitle', 'off');

subplot(2,2,1)
stairs(kT, Xs_k_ot(1,:), 'b-', 'LineWidth', 1.5), grid on
ylabel('x_{s,1}(k)  (°C)'), title('SS Target State — x_{s,1}(k)')
xlim([0, N-1])

subplot(2,2,2)
stairs(kT, Xs_k_ot(2,:), 'b-', 'LineWidth', 1.5), grid on
ylabel('x_{s,2}(k)  (°C)'), title('SS Target State — x_{s,2}(k)')
xlim([0, N-1])

subplot(2,2,3)
stairs(kT, Us_k_ot(1,:), 'r-', 'LineWidth', 1.5), hold on, grid on
yline(100, 'k--', 'LineWidth', 0.8)
yline(  0, 'k--', 'LineWidth', 0.8)
ylabel('u_{s,1}(k)  (%)'), xlabel('Sample k'), title('SS Target Input — u_{s,1}(k)')
ylim([-5, 105]), xlim([0, N-1])

subplot(2,2,4)
stairs(kT, Us_k_ot(2,:), 'r-', 'LineWidth', 1.5), hold on, grid on
yline(100, 'k--', 'LineWidth', 0.8)
yline(  0, 'k--', 'LineWidth', 0.8)
ylabel('u_{s,2}(k)  (%)'), xlabel('Sample k'), title('SS Target Input — u_{s,2}(k)')
ylim([-5, 105]), xlim([0, N-1])

%% ── 16.  Save Experiment Results ─────────────────────────────────────────
save TCL_MPC_OT_ExpResults ...
    Yk_ot Uk_ot Rk_f_ot ek_f_ot Xs_k_ot Us_k_ot ...
    kT Ts Ys Us Xs n_op n_ip N ...
    ISE_T1 ISE_T2 IAE_T1 IAE_T2 ...
    Np Nc Wx Wdelu beta_r alfa_e

fprintf('Results saved to TCL_MPC_OT_ExpResults.mat\n');
fprintf('  Export variables: Yk_ot, Uk_ot, Rk_f_ot, ek_f_ot, Xs_k_ot, Us_k_ot\n');
fprintf('  Timing:           kT  (sample indices 0..N-1),  Ts = %g s\n', Ts);
