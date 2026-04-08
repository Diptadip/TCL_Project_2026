close all; clear all; clc

%% ── 0.  Connect to Arduino ───────────────────────────────────────────────
tclab_N;

%% ── 1.  Load mechanistic model ───────────────────────────────────────────
disp('Loading mechanistic model and computing MPC matrices ...')
load TCL_MechModel_Parameters.mat

A  = TCL.phy;        % Ad  (n_st × n_st)
B  = TCL.gama;       % Bd  (n_st × n_ip)
C  = TCL.C_mat;      % Cd  (n_op × n_st)
D  = TCL.D_mat;      % Dd  (n_op × n_ip)
Ts = TCL.Samp_T;     % 4 s

[n_op, n_st] = size(C);
[~,    n_ip] = size(B);

Ys = TCL.Ys(:);
Us = TCL.Us(:);
Xs = TCL.Xs(:);

fprintf('=== TCL Output-Tracking MPC Experiment (Mechanistic Model) ===\n');
fprintf('  n_st=%d  n_ip=%d  n_op=%d   Ts=%g s\n', n_st, n_ip, n_op, Ts);
fprintf('  Operating point:  T1_ss=%.2f degC   T2_ss=%.2f degC\n', Ys(1), Ys(2));
fprintf('  Operating point:  H1_ss=%.2f %%      H2_ss=%.2f %%\n\n', Us(1), Us(2));

%% ── 2.  MPC Tuning Parameters ────────────────────────────────────────────
Np    = 15;
Nc    = 5;

% Output weights (same as simulation)
Wy    = diag([1,  15 ]);
Wdelu = diag([1.0, 1.0]);

fprintf('=== MPC Tuning ===\n');
fprintf('  Np=%d  Nc=%d\n', Np, Nc);
fprintf('  Wy    = diag([%.1f %.1f])\n',   Wy(1,1),    Wy(2,2));
fprintf('  Wdelu = diag([%.2f %.2f])\n\n', Wdelu(1,1), Wdelu(2,2));

%% ── 3.  Kalman Filter Design ─────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st);
R_kf = TCL.R;
[~, ~, L_kf] = dare(A', C', Q_kf, R_kf);
L_kf = L_kf';
fprintf('Kalman gain L_kf computed (DARE).\n\n');

%% ── 4.  Steady-State Target Computation ──────────────────────────────────
M_ss     = [(eye(n_st) - A), -B ; C, zeros(n_op, n_ip)];
M_ss_inv = pinv(M_ss);

%% ── 5.  Pre-compute QP Matrices (Output Prediction) ─────────────────────
A_aug = [A, B; zeros(n_ip, n_st), eye(n_ip)];
B_aug = [B; eye(n_ip)];
Cy_aug = [C, zeros(n_op, n_ip)];    % output extractor on augmented state

n_aug = n_st + n_ip;
Psi_y   = zeros(Np * n_op, n_aug);
Theta_y = zeros(Np * n_op, Nc * n_ip);

A_aug_pow = eye(n_aug);
for i = 1 : Np
    A_aug_pow       = A_aug_pow * A_aug;
    rows_i          = (i-1)*n_op + 1 : i*n_op;
    Psi_y(rows_i,:) = Cy_aug * A_aug_pow;
    for j = 1 : min(i, Nc)
        cols_j = (j-1)*n_ip + 1 : j*n_ip;
        A_ij   = eye(n_aug);
        for ii = 1 : i-j
            A_ij = A_ij * A_aug;
        end
        Theta_y(rows_i, cols_j) = Cy_aug * A_ij * B_aug;
    end
end

Wy_bar    = kron(eye(Np), Wy);
Wdelu_bar = kron(eye(Nc), Wdelu);

% Hessian H (fixed, pre-computed once)
H_qp = Theta_y' * Wy_bar * Theta_y + Wdelu_bar;
H_qp = (H_qp + H_qp') / 2;         % ensure exact symmetry

% Inequality constraint matrix  M  for cumulative delta-u -> u
M_mat  = kron(tril(ones(Nc)), eye(n_ip));
A_cons = [M_mat; -M_mat];

qp_options = optimoptions('quadprog', 'Display', 'off');

fprintf('QP matrices pre-computed. Ready to start experiment.\n\n');

%% ── 6.  Experiment Timing & Setpoint Schedule ────────────────────────────
N_samp   = 550;
k_warmup = 30;
k_settle = 100;
k_step1  = 101;
k_step2  = 301;

delta_r_step1 = [ +6; -4];
delta_r_step2 = [ -5; +7];

u_H = [100; 100] - Us;       % upper input perturbation bound
u_L = [  0;   0] - Us;       % lower input perturbation bound

fprintf('=== Experiment Setpoints ===\n');
fprintf('  Step 1 (k=%d): T1->%.1f degC  T2->%.1f degC\n', ...
        k_step1, Ys(1)+delta_r_step1(1), Ys(2)+delta_r_step1(2));
fprintf('  Step 2 (k=%d): T1->%.1f degC  T2->%.1f degC\n\n', ...
        k_step2, Ys(1)+delta_r_step2(1), Ys(2)+delta_r_step2(2));
fprintf('  Total duration: %d samples × %g s = %.1f min\n\n', ...
        N_samp, Ts, N_samp*Ts/60);

%% ── 7.  Filter Parameters ────────────────────────────────────────────────
beta_r   = 0.95;   phy_r = beta_r * eye(n_op);
alfa_e   = 0.95;   phy_e = alfa_e * eye(n_op);
rk_f_now = zeros(n_op, 1);   % running filtered setpoint (perturbation)
ek_f_now = zeros(n_op, 1);   % running filtered error

%% ── 8.  MPC and Kalman State Memory ─────────────────────────────────────
xk_hat = zeros(n_st, 1);     % Kalman state estimate (perturbation)
u_prev = zeros(n_ip, 1);     % previous input perturbation (for delta-u)

%% ── 9.  Data Logging ─────────────────────────────────────────────────────
t1s  = [];   t2s  = [];
h1s  = [];   h2s  = [];
R1s  = [];   R2s  = [];
e1s  = [];   e2s  = [];
xs1s = [];   xs2s = [];
us1s = [];   us2s = [];

%% ── 10.  Live Figure Setup ───────────────────────────────────────────────
figure(1)

%% ── 11.  Main Real-Time Loop ─────────────────────────────────────────────
fprintf('Starting MPC-OT experiment.  Total: %.1f min\n', N_samp*Ts/60);
fprintf('Phase 1 (warm-up): heaters OFF for %d samples (%d s)\n\n', ...
        k_warmup, k_warmup*Ts);

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

    % ---- Setpoint pre-filter --------------------------------------------
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
        ek_f_now = zeros(n_op,  1);
        rk_f_now = zeros(n_op,  1);

    else
        % ---- MPC: online constrained QP (quadprog) ----------------------
        x_aug = [xk_hat; u_prev];

        % Output tracking reference over horizon
        R_ref = repmat(rk_f_now, Np, 1);

        % Linear term: f = Theta' * Wy_bar * (Psi*x - Ref)
        f_qp = Theta_y' * Wy_bar * (Psi_y * x_aug - R_ref);

        % Dynamic input bounds (relative to u_prev)
        b_cons = [repmat(u_H - u_prev, Nc, 1);
                  repmat(-u_L + u_prev, Nc, 1)];

        % Solve QP
        [du_full, ~, exitflag] = quadprog(H_qp, f_qp, A_cons, b_cons, ...
                                          [], [], [], [], [], qp_options);

        % Robustness fallback (infeasible / solver fail)
        if isempty(du_full) || exitflag < 0
            du_opt = zeros(n_ip, 1);
        else
            du_opt = du_full(1 : n_ip);  % first move only (receding horizon)
        end

        uk_pert = u_prev + du_opt;
        u_prev  = uk_pert;
    end

    % ---- Absolute heater values -----------------------------------------
    Uk_abs = uk_pert + Us;
    Rk_abs = rk_f_now + Ys;

    % ---- Send to hardware -----------------------------------------------
    h1(Uk_abs(1));
    h2(Uk_abs(2));

    % ---- Kalman time update (prediction step for next sample) -----------
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

    % ---- Phase transition banners ---------------------------------------
    if k == k_warmup + 1
        fprintf('Phase 2: MPC-OT active, settling at Ys ...\n');
    elseif k == k_step1
        fprintf('Phase 3 – Step 1  (T1%+.0f degC, T2%+.0f degC)\n', ...
                delta_r_step1(1), delta_r_step1(2));
    elseif k == k_step2
        fprintf('Phase 4 – Step 2  (T1%+.0f degC, T2%+.0f degC)\n', ...
                delta_r_step2(1), delta_r_step2(2));
    end

    % ---- Live strip chart (every 15 samples ≈ 60 s) ---------------------
    if rem(k, 15) == 0
        n_now = length(t1s);
        kk    = 0 : n_now-1;

        clf(figure(1))
        subplot(2,1,1)
        plot(kk, t1s, 'b-',  kk, R1s, 'b--', ...
             kk, t2s, 'r-',  kk, R2s, 'r--')
        ylabel('Temperature (deg C)')
        legend('T_1','Ref_1','T_2','Ref_2','Location','NorthWest')
        title(sprintf('TCL MPC-OT Experiment  (Mechanistic)  –  k = %d / %d', ...
              k, N_samp))
        grid on

        subplot(2,1,2)
        stairs(kk, h1s, 'b-', 'LineWidth', 2), hold on
        stairs(kk, h2s, 'r-', 'LineWidth', 2)
        yline(  0, 'k--', 'LineWidth', 0.8)
        yline(100, 'k--', 'LineWidth', 0.8)
        ylabel('Heater (%)'), xlabel('Sample k')
        legend('H_1','H_2','Location','NorthWest')
        ylim([-5, 105]), grid on
        drawnow
    end

    % ---- Console output -------------------------------------------------
    elapsed = toc;
    fprintf('k=%3d | T1=%.2f  T2=%.2f | H1=%5.1f%%  H2=%5.1f%% | slack=%.3fs\n', ...
        k, Tk(1), Tk(2), Uk_abs(1), Uk_abs(2), Ts - elapsed);
    pause(max(0.01, Ts - elapsed))

end  % ── end main loop ──────────────────────────────────────────────────

%% ── 12.  Turn heaters off safely ────────────────────────────────────────
h1(0); h2(0);
fprintf('\nExperiment complete. Heaters set to 0%%.\n\n');

%% ── 13.  Build absolute arrays (mirror simulation naming) ───────────────
N   = N_samp;
kT  = (0 : N-1)';
t_sec = kT * Ts;

Yk_ot    = [t1s; t2s];                          % 2 × N  absolute outputs
Uk_ot    = [h1s; h2s];                          % 2 × N  absolute inputs
Rk_f_ot  = [R1s; R2s];                          % 2 × N  absolute reference
ek_f_ot  = [e1s; e2s];                          % 2 × N  filtered error
Xs_k_ot  = [xs1s; xs2s];                        % 2 × N  absolute SS states
Us_k_ot  = [us1s; us2s];                        % 2 × N  absolute SS inputs

%% ── 14.  Performance Metrics ─────────────────────────────────────────────
ISE_T1 = sum((Yk_ot(1,:) - Rk_f_ot(1,:)).^2) * Ts;
ISE_T2 = sum((Yk_ot(2,:) - Rk_f_ot(2,:)).^2) * Ts;
IAE_T1 = sum(abs(Yk_ot(1,:) - Rk_f_ot(1,:))) * Ts;
IAE_T2 = sum(abs(Yk_ot(2,:) - Rk_f_ot(2,:))) * Ts;

fprintf('=== Tracking Performance ===\n');
fprintf('  ISE T1 = %8.2f   |   ISE T2 = %8.2f\n', ISE_T1, ISE_T2);
fprintf('  IAE T1 = %8.2f   |   IAE T2 = %8.2f\n\n', IAE_T1, IAE_T2);

%% ── 15.  Post-Experiment Plots ───────────────────────────────────────────

% ── FIGURE 2: Output profiles ────────────────────────────────────────────
figure(2), clf
subplot(2,1,1)
plot(kT, Yk_ot(1,:), 'b-', 'LineWidth', 1.5), hold on
plot(kT, Rk_f_ot(1,:), 'b--', 'LineWidth', 1.0)
grid on, ylabel('T_1  (deg C)')
title('Output-Tracking MPC (Mech) – Experiment – Output Profiles')
legend('Y_1(k)', 'R_1(k)', 'Location', 'Best')
xline(k_step1-1, 'k:', 'Step 1', 'LabelVerticalAlignment', 'bottom')
xline(k_step2-1, 'k:', 'Step 2', 'LabelVerticalAlignment', 'bottom')

subplot(2,1,2)
plot(kT, Yk_ot(2,:), 'r-', 'LineWidth', 1.5), hold on
plot(kT, Rk_f_ot(2,:), 'r--', 'LineWidth', 1.0)
grid on, xlabel('Sample k'), ylabel('T_2  (deg C)')
legend('Y_2(k)', 'R_2(k)', 'Location', 'Best')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

% ── FIGURE 3: Heater inputs ───────────────────────────────────────────────
figure(3), clf
subplot(2,1,1)
stairs(kT, Uk_ot(1,:)', 'b-', 'LineWidth', 1.5), hold on
yline(  0, 'k--', 'LineWidth', 0.8), yline(100, 'k--', 'LineWidth', 0.8)
grid on, ylabel('U_1(k)  (% heater)')
title('Output-Tracking MPC (Mech) – Experiment – Heater Inputs')
legend('U_1(k)', 'Bounds', 'Location', 'Best'), ylim([-5, 105])
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

subplot(2,1,2)
stairs(kT, Uk_ot(2,:)', 'r-', 'LineWidth', 1.5), hold on
yline(  0, 'k--', 'LineWidth', 0.8), yline(100, 'k--', 'LineWidth', 0.8)
grid on, xlabel('Sample k'), ylabel('U_2(k)  (% heater)')
legend('U_2(k)', 'Bounds', 'Location', 'Best'), ylim([-5, 105])
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

% ── FIGURE 4: Filtered tracking error ────────────────────────────────────
figure(4), clf
subplot(2,1,1)
stairs(kT, ek_f_ot(1,:)', 'b-', 'LineWidth', 1.5), hold on
yline(0, 'k--', 'LineWidth', 0.8), grid on
ylabel('e_{f,1}(k)  (deg C)')
title('Output-Tracking MPC (Mech) – Experiment – Filtered Error e_f(k)')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

subplot(2,1,2)
stairs(kT, ek_f_ot(2,:)', 'r-', 'LineWidth', 1.5), hold on
yline(0, 'k--', 'LineWidth', 0.8), grid on
xlabel('Sample k'), ylabel('e_{f,2}(k)  (deg C)')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

% ── FIGURE 5: Steady-state targets ───────────────────────────────────────
figure(5), clf
subplot(2,1,1)
plot(kT, Xs_k_ot(1,:), 'b-', 'LineWidth', 1.5), hold on
plot(kT, Xs_k_ot(2,:), 'r-', 'LineWidth', 1.5)
grid on, ylabel('x_s(k)  (deg C)')
title('Output-Tracking MPC (Mech) – Experiment – SS Target States x_s(k)')
legend('x_{s,1}(k)', 'x_{s,2}(k)', 'Location', 'Best')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

subplot(2,1,2)
stairs(kT, Us_k_ot(1,:)', 'b-', 'LineWidth', 1.5), hold on
stairs(kT, Us_k_ot(2,:)', 'r-', 'LineWidth', 1.5)
grid on, xlabel('Sample k'), ylabel('u_s(k)  (% heater)')
title('Output-Tracking MPC (Mech) – Experiment – SS Target Inputs u_s(k)')
legend('u_{s,1}(k)', 'u_{s,2}(k)', 'Location', 'Best')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

%% ── 16.  Save results ────────────────────────────────────────────────────
save TCL_MPC_OT_ExpResults ...
    kT t_sec Yk_ot Uk_ot Rk_f_ot ek_f_ot Xs_k_ot Us_k_ot ...
    ISE_T1 ISE_T2 IAE_T1 IAE_T2 ...
    Ts Ys Us Xs n_op n_ip N ...
    Np Nc Wy Wdelu beta_r alfa_e

fprintf('Saved: TCL_MPC_OT_ExpResults.mat\n');
fprintf('All done.\n');