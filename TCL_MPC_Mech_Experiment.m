% =========================================================================
close all; clear all; clc

global TCL

%% ── 0.  Connect to Arduino ───────────────────────────────────────────────
tclab_N;

%% ── 1.  Load mechanistic model parameters ────────────────────────────────
disp('Loading mechanistic model and computing MPC parameters ...')
load TCL_MechModel_Parameters.mat    % TCL struct

A  = TCL.phy;       % [n_st × n_st]
B  = TCL.gama;      % [n_st × n_ip]
C  = TCL.C_mat;     % [n_op × n_st]
D  = TCL.D_mat;     % [n_op × n_ip]
Ts = TCL.Samp_T;    % 4 s

[n_op, n_st] = size(C);
[~,    n_ip] = size(B);

Ys = TCL.Ys(:);    % steady-state outputs  [deg C]
Us = TCL.Us(:);    % steady-state inputs   [%]
Xs = TCL.Xs(:);    % steady-state states   [deg C]

fprintf('=== TCL State-Tracking MPC (Mechanistic) – EXPERIMENT ===\n');
fprintf('  State dim n_st = %d,  Inputs n_ip = %d,  Outputs n_op = %d\n', ...
        n_st, n_ip, n_op);
fprintf('  Operating point:  T1_ss = %.2f degC   T2_ss = %.2f degC\n', Ys(1), Ys(2));
fprintf('  Operating point:  H1_ss = %.2f %%      H2_ss = %.2f %%\n\n', Us(1), Us(2));

%% ── 2.  MPC Tuning Parameters ────────────────────────────────────────────
Np    = 15;
Nc    = 5;

Wx    = diag([15,  1 ]);              % state tracking weight (n_st × n_st)
Wdelu = diag([1.0, 1.0]);             % input-move weight     (n_ip × n_ip)

fprintf('=== MPC Tuning ===\n');
fprintf('  Np = %d,  Nc = %d,  Ts = %g s\n', Np, Nc, Ts);
fprintf('  Wx    = diag([%.1f %.1f])\n', Wx(1,1), Wx(2,2));
fprintf('  Wdelu = diag([%.1f %.1f])\n\n', Wdelu(1,1), Wdelu(2,2));

%% ── 3.  Kalman Filter Design ─────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st);
R_kf = TCL.R;
[~, ~, L_kf] = dare(A', C', Q_kf, R_kf);
L_kf = L_kf';                         % [n_st × n_op]

fprintf('Kalman gain L_kf computed (DARE).\n\n');

%% ── 4.  Steady-State Target Gain ─────────────────────────────────────────
M_ss     = [(eye(n_st) - A), -B; C, zeros(n_op, n_ip)];
M_ss_inv = pinv(M_ss);

fprintf('Steady-state target gain M_ss_inv pre-computed.\n\n');

%% ── 5.  Pre-compute QP Matrices (State Prediction) ───────────────────────
%  Augmented state: xi(k) = [x_hat(k); u(k-1)]
A_aug = [A, B; zeros(n_ip, n_st), eye(n_ip)];
B_aug = [B; eye(n_ip)];
Cx_aug = [eye(n_st), zeros(n_st, n_ip)];   % extracts state from xi

n_aug = n_st + n_ip;

Psi_x   = zeros(Np*n_st, n_aug);
Theta_x = zeros(Np*n_st, Nc*n_ip);

A_aug_pow = eye(n_aug);
for i = 1:Np
    A_aug_pow         = A_aug_pow * A_aug;
    rows_i            = (i-1)*n_st + 1 : i*n_st;
    Psi_x(rows_i, :)  = Cx_aug * A_aug_pow;

    for j = 1:min(i, Nc)
        cols_j   = (j-1)*n_ip + 1 : j*n_ip;
        A_aug_ij = eye(n_aug);
        for ii = 1:i-j
            A_aug_ij = A_aug_ij * A_aug;
        end
        Theta_x(rows_i, cols_j) = Cx_aug * A_aug_ij * B_aug;
    end
end

Wx_bar    = kron(eye(Np), Wx);
Wdelu_bar = kron(eye(Nc), Wdelu);

% Hessian for quadprog (fixed across all steps)
H_qp = Theta_x' * Wx_bar * Theta_x + Wdelu_bar;
H_qp = (H_qp + H_qp') / 2;           % enforce exact symmetry

% Constraint matrix M for cumulative input increments
M_mat  = kron(tril(ones(Nc)), eye(n_ip));
A_cons = [M_mat; -M_mat];

qp_options = optimoptions('quadprog', 'Display', 'off');

fprintf('QP matrices pre-computed.  Ready to start experiment.\n\n');

%% ── 6.  Experiment Timing ────────────────────────────────────────────────
N_samp   = 550;
k_warmup = 30;
k_step1  = 101;
k_step2  = 301;

delta_r_step1 = [ +6; -4 ];
delta_r_step2 = [ -5; +7 ];

u_H = [100; 100] - Us;   % upper input-perturbation bound
u_L = [  0;   0] - Us;   % lower input-perturbation bound

fprintf('=== Experiment Setpoints ===\n');
fprintf('  Step 1 (k=%d): T1 -> %.1f degC,  T2 -> %.1f degC\n', ...
        k_step1, Ys(1)+delta_r_step1(1), Ys(2)+delta_r_step1(2));
fprintf('  Step 2 (k=%d): T1 -> %.1f degC,  T2 -> %.1f degC\n\n', ...
        k_step2, Ys(1)+delta_r_step2(1), Ys(2)+delta_r_step2(2));
fprintf('  Total duration: %d samples × %g s = %.1f min\n\n', ...
        N_samp, Ts, N_samp*Ts/60);

%% ── 7.  Filter Parameters ────────────────────────────────────────────────
beta_r   = 0.95;   phy_r = beta_r * eye(n_op);
alfa_e   = 0.95;   phy_e = alfa_e * eye(n_op);

rk_f_now = zeros(n_op, 1);
ek_f_now = zeros(n_op, 1);

%% ── 8.  MPC and Kalman State Memory ─────────────────────────────────────
xk_hat = zeros(n_st, 1);   % Kalman state estimate (perturbation)
u_prev = zeros(n_ip, 1);   % previous input perturbation

%% ── 9.  Data Logging Arrays ──────────────────────────────────────────────
t1s  = [];  t2s  = [];
h1s  = [];  h2s  = [];
R1s  = [];  R2s  = [];
e1s  = [];  e2s  = [];
xs1s = [];  xs2s = [];
us1s = [];  us2s = [];

%% ── 10.  Live Figure ─────────────────────────────────────────────────────
figure(1); clf;

%% ── 11.  Main Real-Time Loop ─────────────────────────────────────────────
fprintf('Starting experiment.  Phase 1 (warm-up): heaters OFF for %d s\n\n', ...
        k_warmup*Ts);

for k = 1:N_samp
    tic

    % ---- Read hardware --------------------------------------------------
    Tk = [T1C(); T2C()];          % absolute temperatures [deg C]
    yk_meas = Tk - Ys;            % output perturbation

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

    % ---- Innovation filter (for logging ek_f) ---------------------------
    raw_ek   = rk_f_now - yk_meas;
    ek_f_now = phy_e * ek_f_now + (eye(n_op) - phy_e) * raw_ek;

    % ---- Kalman measurement correction ----------------------------------
    y_pred = C * xk_hat;
    inno   = yk_meas - y_pred;
    xk_hat = xk_hat + L_kf * inno;   % corrected state estimate

    % ---- Steady-state targets -------------------------------------------
    rhs  = [zeros(n_st, 1); rk_f_now];
    xus  = M_ss_inv * rhs;
    xs_t = xus(1:n_st);
    us_t = xus(n_st+1:end);

    % ---- Warm-up phase: heaters OFF -------------------------------------
    if k <= k_warmup
        uk_pert  = zeros(n_ip, 1);
        u_prev   = zeros(n_ip, 1);
        xk_hat   = zeros(n_st, 1);
        ek_f_now = zeros(n_op, 1);
        rk_f_now = zeros(n_op, 1);
    else
        % ---- MPC constrained QP -----------------------------------------
        x_aug = [xk_hat; u_prev];

        % State-tracking reference (tiled over prediction horizon)
        X_ref = repmat(xs_t, Np, 1);

        % Linear QP cost term
        f_qp = Theta_x' * Wx_bar * (Psi_x * x_aug - X_ref);

        % Dynamic input-constraint bounds (relative to u_prev)
        b_cons = [repmat(u_H - u_prev, Nc, 1);
                  repmat(-u_L + u_prev, Nc, 1)];

        % Solve QP
        [du_full, ~, exitflag] = quadprog(H_qp, f_qp, A_cons, b_cons, ...
                                           [], [], [], [], [], qp_options);

        % Robustness fallback
        if isempty(du_full) || exitflag < 0
            du_opt = zeros(n_ip, 1);
        else
            du_opt = du_full(1:n_ip);   % first move only (receding horizon)
        end

        uk_pert = u_prev + du_opt;

        % Hard clamp (anti-windup safety net)
        uk_pert = min(max(uk_pert, u_L), u_H);

        u_prev = uk_pert;
    end

    % ---- Absolute heater commands ----------------------------------------
    Uk_abs = uk_pert + Us;
    Rk_abs = rk_f_now + Ys;

    % ---- Send to hardware -----------------------------------------------
    h1(Uk_abs(1));
    h2(Uk_abs(2));

    % ---- Kalman time-update (prediction for next step) ------------------
    xk_hat = A * xk_hat + B * uk_pert;

    % ---- Log data -------------------------------------------------------
    t1s  = [t1s,  Tk(1)];
    t2s  = [t2s,  Tk(2)];
    h1s  = [h1s,  Uk_abs(1)];
    h2s  = [h2s,  Uk_abs(2)];
    R1s  = [R1s,  Rk_abs(1)];
    R2s  = [R2s,  Rk_abs(2)];
    e1s  = [e1s,  ek_f_now(1)];
    e2s  = [e2s,  ek_f_now(2)];
    xs1s = [xs1s, xs_t(1) + Xs(1)];
    xs2s = [xs2s, xs_t(2) + Xs(2)];
    us1s = [us1s, us_t(1) + Us(1)];
    us2s = [us2s, us_t(2) + Us(2)];

    % ---- Phase transition banners ----------------------------------------
    if k == k_warmup + 1
        fprintf('Phase 2: MPC active, settling at operating point ...\n');
    elseif k == k_step1
        fprintf('Phase 3: Step 1  (T1 %+.0f degC, T2 %+.0f degC)\n', ...
                delta_r_step1(1), delta_r_step1(2));
    elseif k == k_step2
        fprintf('Phase 4: Step 2  (T1 %+.0f degC, T2 %+.0f degC)\n', ...
                delta_r_step2(1), delta_r_step2(2));
    end

    % ---- Live strip chart (every 15 samples ≈ 60 s) ---------------------
    if rem(k, 15) == 0
        n_now = length(t1s);
        kk    = 0:n_now-1;

        clf(figure(1))
        subplot(2,1,1)
        plot(kk, t1s, 'b-',  kk, R1s, 'b--', ...
             kk, t2s, 'r-',  kk, R2s, 'r--', 'LineWidth', 1.5)
        ylabel('Temperature (deg C)')
        legend('T_1','Ref_1','T_2','Ref_2','Location','NorthWest')
        title(sprintf('TCL MPC-Mech (State-Tracking) Experiment  –  k = %d / %d', ...
                      k, N_samp))
        grid on
        xline(k_step1-1, 'k:', 'Step 1', 'LabelVerticalAlignment','bottom')
        xline(k_step2-1, 'k:', 'Step 2', 'LabelVerticalAlignment','bottom')

        subplot(2,1,2)
        stairs(kk, h1s, 'b-', 'LineWidth', 1.5), hold on
        stairs(kk, h2s, 'r-', 'LineWidth', 1.5)
        yline(5,  'k--', 'LineWidth', 0.8)
        yline(80, 'k--', 'LineWidth', 0.8)
        ylabel('Heater (%)'), xlabel('Sample k')
        legend('H_1','H_2','Safety bounds','Location','NorthWest')
        grid on
        xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')
        drawnow
    end

    % ---- Console output -------------------------------------------------
    elapsed = toc;
    fprintf('k=%3d | T1=%.2f T2=%.2f | H1=%.1f%% H2=%.1f%% | slack=%.3fs\n', ...
            k, Tk(1), Tk(2), Uk_abs(1), Uk_abs(2), Ts - elapsed);

    pause(max(0.01, Ts - elapsed))

end  % main loop

%% ── 12.  Turn Heaters Off ────────────────────────────────────────────────
h1(0); h2(0);
fprintf('\nExperiment complete. Heaters set to 0%%.\n\n');

%% ── 13.  Build Result Arrays ─────────────────────────────────────────────
kT    = (0:N_samp-1)';
t_sec = kT * Ts;

Yk   = [t1s; t2s];
Uk   = [h1s; h2s];
Rk_f = [R1s; R2s];
ek_f = [e1s; e2s];
Xs_k = [xs1s; xs2s];
Us_k = [us1s; us2s];

%% ── 14.  Post-Experiment Plots ───────────────────────────────────────────
% -- Figure 2: Temperature profiles --
figure(2); clf
subplot(211)
plot(kT, Yk(1,:), 'b-', 'LineWidth', 1.5), hold on
plot(kT, Rk_f(1,:), 'b--', 'LineWidth', 1.0)
grid on, ylabel('T_1  (deg C)')
title('State-Tracking MPC – Experiment – Output Profiles')
legend('Y_1(k)', 'R_1(k)', 'Location','Best')
xline(k_step1-1, 'k:', 'Step 1', 'LabelVerticalAlignment','bottom')
xline(k_step2-1, 'k:', 'Step 2', 'LabelVerticalAlignment','bottom')

subplot(212)
plot(kT, Yk(2,:), 'r-', 'LineWidth', 1.5), hold on
plot(kT, Rk_f(2,:), 'r--', 'LineWidth', 1.0)
grid on, xlabel('Sample k'), ylabel('T_2  (deg C)')
legend('Y_2(k)', 'R_2(k)', 'Location','Best')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

% -- Figure 3: Heater inputs --
figure(3); clf
subplot(211)
stairs(kT, Uk(1,:)', 'b-', 'LineWidth', 1.5), hold on
yline(5,  'k--', 'LineWidth', 0.8), yline(80, 'k--', 'LineWidth', 0.8)
grid on, ylabel('U_1(k)  (% heater)')
title('State-Tracking MPC – Experiment – Heater Inputs')
legend('U_1(k)', 'Safety bounds', 'Location','Best')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

subplot(212)
stairs(kT, Uk(2,:)', 'r-', 'LineWidth', 1.5), hold on
yline(5,  'k--', 'LineWidth', 0.8), yline(80, 'k--', 'LineWidth', 0.8)
grid on, xlabel('Sample k'), ylabel('U_2(k)  (% heater)')
legend('U_2(k)', 'Safety bounds', 'Location','Best')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

% -- Figure 4: Filtered error --
figure(4); clf
subplot(211)
stairs(kT, ek_f(1,:)', 'b-', 'LineWidth', 1.5), hold on
yline(0, 'k--', 'LineWidth', 0.8), grid on
ylabel('e_{f,1}(k)  (deg C)')
title('State-Tracking MPC – Experiment – Filtered Error e_f(k)')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

subplot(212)
stairs(kT, ek_f(2,:)', 'r-', 'LineWidth', 1.5), hold on
yline(0, 'k--', 'LineWidth', 0.8), grid on
xlabel('Sample k'), ylabel('e_{f,2}(k)  (deg C)')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

% -- Figure 5: Steady-state targets --
figure(5); clf
subplot(211)
plot(kT, Xs_k(1,:), 'b-', 'LineWidth', 1.5), hold on
plot(kT, Xs_k(2,:), 'r-', 'LineWidth', 1.5)
grid on, ylabel('x_s(k)  (deg C)')
title('State-Tracking MPC – Experiment – SS Target States x_s(k)')
legend('x_{s,1}(k)', 'x_{s,2}(k)', 'Location','Best')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

subplot(212)
stairs(kT, Us_k(1,:)', 'b-', 'LineWidth', 1.5), hold on
stairs(kT, Us_k(2,:)', 'r-', 'LineWidth', 1.5)
grid on, xlabel('Sample k'), ylabel('u_s(k)  (% heater)')
title('State-Tracking MPC – Experiment – SS Target Inputs u_s(k)')
legend('u_{s,1}(k)', 'u_{s,2}(k)', 'Location','Best')
xline(k_step1-1, 'k:'), xline(k_step2-1, 'k:')

%% ── 15.  Save Results ────────────────────────────────────────────────────
save TCL_MPC_Mech_ExpResults  kT t_sec Yk Uk Rk_f ek_f Xs_k Us_k
fprintf('Results saved to TCL_MPC_Mech_ExpResults.mat\n');
fprintf('\nAll done.\n');
