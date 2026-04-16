close all; clear all; clc

%% ── 0.  Connect to Arduino ───────────────────────────────────────────────
tclab_N;

%% ── 1.  Load models ──────────────────────────────────────────────────────
disp('Loading model and MPC parameters ...')
load TCL_MechModel_Parameters.mat
load TCL_BlackBox_OE_SS.mat

A  = idmod.phy;
B  = idmod.gama;
C  = idmod.C_mat;
D  = idmod.D_mat;
Ts = idmod.samp_T;   % 4 s

[n_op, n_st] = size(C);
[~,    n_ip] = size(B);

Ys = TCL.Ys(:);
Us = TCL.Us(:);

fprintf('Operating point:  T1_ss=%.2f degC   T2_ss=%.2f degC\n', Ys(1), Ys(2));
fprintf('                  H1_ss=%.2f %%       H2_ss=%.2f %%\n\n', Us(1), Us(2));

%% ── 2.  MPC Tuning Parameters ────────────────────────────────────────────
Np    = 15;
Nc    = 5;
Wx    = diag([15.0, 15.0]);
Wdelu = diag([1.0, 1.0]);

%% ── 3.  Kalman Filter ────────────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st);
R_kf = TCL.R;
[~, ~, L_kf] = dare(A', C', Q_kf, R_kf);
L_kf = L_kf';

%% ── 4.  Steady-State Target Gain ─────────────────────────────────────────
M_ss     = [(eye(n_st)-A), -B; C, zeros(n_op, n_ip)];
M_ss_inv = pinv(M_ss);

%% ── 5.  Pre-compute MPC Gain ─────────────────────────────────────────────
A_aug = [A, B; zeros(n_ip, n_st), eye(n_ip)];
B_aug = [B; eye(n_ip)];
C_aug = [C, zeros(n_op, n_ip)];
n_aug = n_st + n_ip;

Psi   = zeros(Np*n_op, n_aug);
Theta = zeros(Np*n_op, Nc*n_ip);

A_aug_pow = eye(n_aug);
for i = 1:Np
    A_aug_pow      = A_aug_pow * A_aug;
    rows_i         = (i-1)*n_op + 1 : i*n_op;
    Psi(rows_i,:)  = C_aug * A_aug_pow;
    for j = 1:min(i, Nc)
        cols_j = (j-1)*n_ip + 1 : j*n_ip;
        A_ij   = eye(n_aug);
        for ii = 1:i-j
            A_ij = A_ij * A_aug;
        end
        Theta(rows_i, cols_j) = C_aug * A_ij * B_aug;
    end
end

Wx_bar    = kron(eye(Np), Wx);
Wdelu_bar = kron(eye(Nc), Wdelu);
H_qp      = Theta' * Wx_bar * Theta + Wdelu_bar;
K_mpc     = (H_qp \ (Theta' * Wx_bar));
K_mpc     = K_mpc(1:n_ip, :);   % first move block only

disp('MPC gain pre-computed. Ready to start experiment.')

%% ── 6.  Experiment Timing ────────────────────────────────────────────────
N_samp   = 550;
k_warmup = 30;
k_settle = 100;
k_step1  = 101;
k_step2  = 301;

delta_r_step1 = [ +4; -2 ];
delta_r_step2 = [ -2; +4 ];

u_H = [100; 100] - Us;
u_L = [  0;   0] - Us;

fprintf('=== Experiment Setpoints ===\n');
fprintf('  Step 1: T1->%.1f degC  T2->%.1f degC\n', Ys(1)+delta_r_step1(1), Ys(2)+delta_r_step1(2));
fprintf('  Step 2: T1->%.1f degC  T2->%.1f degC\n\n', Ys(1)+delta_r_step2(1), Ys(2)+delta_r_step2(2));

%% ── 7.  Filter Parameters (SAME as PI servo experiment) ──────────────────
beta_r   = 0.95;   phy_r = beta_r * eye(n_op);
alfa_e   = 0.95;   phy_e = alfa_e * eye(n_op);
rk_f_now = zeros(n_op, 1);
ek_f_now = zeros(n_op, 1);

%% ── 8.  MPC and Kalman state memory ─────────────────────────────────────
xk_hat = zeros(n_st, 1);   % Kalman state estimate
u_prev = zeros(n_ip, 1);   % previous input perturbation

%% ── 9.  Data logging ────────────────────────────────────────────────────
t1s  = [];  t2s  = [];
h1s  = [];  h2s  = [];
R1s  = [];  R2s  = [];
e1s  = [];  e2s  = [];
xs1s = [];  xs2s = [];
us1s = [];  us2s = [];

%% ── 10.  Live figure ─────────────────────────────────────────────────────
figure(1)

%% ── 11.  Main real-time loop ─────────────────────────────────────────────
fprintf('Starting MPC experiment. Total: %.1f min\n', N_samp*Ts/60);
fprintf('Phase 1 (warm-up): heaters OFF for %d s\n\n', k_warmup*Ts);

for k = 1:N_samp
    tic

    % ---- Read hardware --------------------------------------------------
    Tk    = zeros(n_op,1);
    Tk(1) = T1C();
    Tk(2) = T2C();
    yk_meas = Tk - Ys;    % output perturbation

    % ---- Setpoint schedule ----------------------------------------------
    if k < k_step1
        setpt = zeros(n_op,1);
    elseif k < k_step2
        setpt = delta_r_step1;
    else
        setpt = delta_r_step2;
    end

    % ---- Setpoint filter ------------------------------------------------
    rk_f_now = phy_r * rk_f_now + (eye(n_op) - phy_r) * setpt;

    % ---- Innovation filter (for logging ef) ----------------------------
    raw_ek   = rk_f_now - yk_meas;
    ek_f_now = phy_e * ek_f_now + (eye(n_op) - phy_e) * raw_ek;

    % ---- Kalman correction ----------------------------------------------
    y_pred = C * xk_hat;
    inno   = yk_meas - y_pred;
    xk_hat = xk_hat + L_kf * inno;   % corrected estimate

    % ---- Steady-state targets -------------------------------------------
    rhs  = [zeros(n_st,1); rk_f_now];
    xus  = M_ss_inv * rhs;
    xs_t = xus(1:n_st);
    us_t = xus(n_st+1:end);

    % ---- Warm-up: heaters OFF -------------------------------------------
    if k <= k_warmup
        uk_pert  = zeros(n_ip,1);
        u_prev   = zeros(n_ip,1);
        xk_hat   = zeros(n_st,1);
        ek_f_now = zeros(n_op,1);
        rk_f_now = zeros(n_op,1);
    else
        % ---- MPC law ---------------------------------------------------
        x_aug  = [xk_hat; u_prev];
        R_vec  = repmat(rk_f_now, Np, 1);
        du_opt = K_mpc * (R_vec - Psi * x_aug);
        u_new  = u_prev + du_opt;

        % ---- Anti-windup -----------------------------------------------
        u_new = min(max(u_new, u_L), u_H);

        uk_pert = u_new;
        u_prev  = u_new;
    end

    % ---- Absolute heater values ----------------------------------------
    Uk_abs = uk_pert + Us;
    Rk_abs = rk_f_now + Ys;

    % ---- Send to hardware -----------------------------------------------
    h1(Uk_abs(1));
    h2(Uk_abs(2));

    % ---- Kalman prediction for next step --------------------------------
    xk_hat = A * xk_hat + B * uk_pert;   % time update

    % ---- Log -----------------------------------------------------------
    t1s  = [t1s,  Tk(1)];
    t2s  = [t2s,  Tk(2)];
    h1s  = [h1s,  Uk_abs(1)];
    h2s  = [h2s,  Uk_abs(2)];
    R1s  = [R1s,  Rk_abs(1)];
    R2s  = [R2s,  Rk_abs(2)];
    e1s  = [e1s,  ek_f_now(1)];
    e2s  = [e2s,  ek_f_now(2)];
    xs1s = [xs1s, xs_t(1)];
    xs2s = [xs2s, xs_t(2)];
    us1s = [us1s, us_t(1)];
    us2s = [us2s, us_t(2)];

    % ---- Phase banners -------------------------------------------------
    if k == k_warmup + 1
        fprintf('Phase 2: MPC active, settling at Ys ...\n');
    elseif k == k_step1
        fprintf('Phase 3: Step 1  (T1%+.0f degC, T2%+.0f degC)\n', delta_r_step1(1), delta_r_step1(2));
    elseif k == k_step2
        fprintf('Phase 4: Step 2  (T1%+.0f degC, T2%+.0f degC)\n', delta_r_step2(1), delta_r_step2(2));
    end

    % ---- Live strip chart (every 15 samples = 60 s) --------------------
    if rem(k, 15) == 0
        n_now = length(t1s);
        kk    = 0:n_now-1;

        clf(figure(1))
        subplot(2,1,1)
        plot(kk, t1s, 'r-',  kk, R1s, 'r--', ...
             kk, t2s, 'b-',  kk, R2s, 'b--')
        ylabel('Temperature (deg C)')
        legend('T_1','Ref_1','T_2','Ref_2','Location','NorthWest')
        title(sprintf('TCL MPC-BB Experiment  –  k = %d / %d', k, N_samp))
        grid on

        subplot(2,1,2)
        stairs(kk, h1s, 'r-', 'LineWidth', 2), hold on
        stairs(kk, h2s, 'b-', 'LineWidth', 2)
        yline(5,  'k--', 'LineWidth', 0.8)
        yline(80, 'k--', 'LineWidth', 0.8)
        ylabel('Heater (%)')
        xlabel('Sample k')
        legend('H_1','H_2','Location','NorthWest')
        grid on
        drawnow
    end

    % ---- Console output -------------------------------------------------
    elapsed = toc;
    fprintf('k=%3d | T1=%.2f T2=%.2f | H1=%.1f%% H2=%.1f%% | slack=%.3fs\n', ...
        k, Tk(1), Tk(2), Uk_abs(1), Uk_abs(2), Ts - elapsed);
    pause(max(0.01, Ts - elapsed))

end  % main loop

%% ── 12.  Turn heaters off ────────────────────────────────────────────────
disp('Turning heaters off ...')
h1(0);
h2(0);
disp('MPC experiment complete.')

%% ── 13.  Post-experiment time axis ──────────────────────────────────────
n_log = length(t1s);
kT    = 0:n_log-1;
tf    = (N_samp-1);   % last sample index

%% ── 14.  Figure 2 – Output profiles Yi(k) vs k and Ri(k) vs k ──────────
figure(2), clf
subplot(211)
plot(kT, t1s, 'b-', kT, R1s, 'b--', 'LineWidth',1.5), grid on
ylabel('T_1  (deg C)')
legend('Y_1(k)','R_1(k)','Location','NorthWest')
title('MPC BB Experiment – Temperature Profiles')
xline(k_step1-1,'k:','Step 1','LabelVerticalAlignment','bottom')
xline(k_step2-1,'k:','Step 2','LabelVerticalAlignment','bottom')

subplot(212)
plot(kT, t2s, 'r-', kT, R2s, 'r--', 'LineWidth',1.5), grid on
ylabel('T_2  (deg C)')
xlabel('Sample k')
legend('Y_2(k)','R_2(k)','Location','NorthWest')
xline(k_step1-1,'k:')
xline(k_step2-1,'k:')

%% ── 15.  Figure 3 – Heater inputs Ui(k) vs k ────────────────────────────
figure(3), clf
subplot(211)
stairs(kT, h1s, 'b-', 'LineWidth',2), hold on
yline(5,  'k--', 'LineWidth',0.8)
yline(80, 'k--', 'LineWidth',0.8)
ylabel('U_1(k)  (%)')
title('MPC BB Experiment – Heater Inputs')
legend('U_1(k)','Safety bounds','Location','NorthWest')
grid on
xline(k_step1-1,'k:')
xline(k_step2-1,'k:')

subplot(212)
stairs(kT, h2s, 'r-', 'LineWidth',2), hold on
yline(5,  'k--', 'LineWidth',0.8)
yline(80, 'k--', 'LineWidth',0.8)
ylabel('U_2(k)  (%)')
xlabel('Sample k')
legend('U_2(k)','Safety bounds','Location','NorthWest')
grid on
xline(k_step1-1,'k:')
xline(k_step2-1,'k:')

%% ── 16.  Figure 4 – Filtered error ef(k) vs k ───────────────────────────
figure(4), clf
subplot(211)
stairs(kT, e1s, 'b-', 'LineWidth',1.5), hold on
yline(0,'k--','LineWidth',0.8), grid on
ylabel('e_{f,1}(k)  (deg C)')
title('MPC BB Experiment – Filtered Error e_f(k)')
xline(k_step1-1,'k:')
xline(k_step2-1,'k:')

subplot(212)
stairs(kT, e2s, 'r-', 'LineWidth',1.5), hold on
yline(0,'k--','LineWidth',0.8), grid on
ylabel('e_{f,2}(k)  (deg C)')
xlabel('Sample k')
xline(k_step1-1,'k:')
xline(k_step2-1,'k:')

%% ── 17.  Figure 5 – SS targets xs(k) and us(k) vs k ────────────────────
figure(5), clf
subplot(211)
plot(kT, xs1s, 'b-', kT, xs2s, 'r-', 'LineWidth',1.5), grid on
ylabel('x_s(k)  (perturbation)')
title('MPC BB Experiment – SS Target States x_s(k)')
legend('x_{s,1}(k)','x_{s,2}(k)','Location','Best')
xline(k_step1-1,'k:')
xline(k_step2-1,'k:')

subplot(212)
stairs(kT, us1s+Us(1), 'b-', 'LineWidth',1.5), hold on
stairs(kT, us2s+Us(2), 'r-', 'LineWidth',1.5)
grid on
ylabel('u_s(k)  (% heater)')
xlabel('Sample k')
title('SS Target Inputs u_s(k)')
legend('u_{s,1}(k)','u_{s,2}(k)','Location','Best')
xline(k_step1-1,'k:')
xline(k_step2-1,'k:')

%% ── 18.  Save results ────────────────────────────────────────────────────
save TCL_MPC_BB_ExpResults  kT t1s t2s h1s h2s R1s R2s e1s e2s xs1s xs2s us1s us2s
fprintf('\nResults saved to TCL_MPC_BB_ExpResults.mat\n');
