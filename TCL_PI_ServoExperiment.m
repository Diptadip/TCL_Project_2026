close all; clear all; clc

%% ── 0.  Connect to Arduino ───────────────────────────────────────────────
tclab_N;

%% ── 1.  Load PI parameters ───────────────────────────────────────────────
disp('Loading PI controller parameters ...')
load TCL_MechModel_Parameters.mat

if ~isfield(TCL, 'PI')
    error('TCL.PI not found. Run TCL_PI_Controller_Design.m first.');
end

Ys = TCL.Ys;       % steady-state temperatures (deg C), 2×1
Us = TCL.Us;       % steady-state heater values (%),    2×1
Ts = TCL.Samp_T;   % 4 s

n_ip = 2;
n_op = 2;

fprintf('Operating point:  T1_ss = %.2f deg C   T2_ss = %.2f deg C\n',  Ys(1), Ys(2));
fprintf('                  H1_ss = %.2f %%        H2_ss = %.2f %%\n\n', Us(1), Us(2));

%% ── 2.  Experiment timing ────────────────────────────────────────────────
N_samp   = 550;     % 550 × 4 s = 2200 s ≈ 36.7 min

k_warmup = 30;      % k = 1  … 30  : heaters OFF   (120 s)
k_settle = 100;     % k = 31 … 100 : PI on, rk = 0  (280 s)
k_step1  = 101;     % k = 101… 300 : Step 1         (800 s)
k_step2  = 301;     % k = 301… 550 : Step 2         (1000 s)

%% ── 3.  Setpoint perturbations (validated by nonlinear simulation) ───────
delta_r_step1 = [ +6 ; -4 ];   % T1 UP +6 deg C,  T2 DOWN -4 deg C
delta_r_step2 = [ -5 ; +7 ];   % T1 DOWN -5 deg C, T2 UP  +7 deg C

fprintf('=== Servo Experiment Setpoints ===\n');
fprintf('  Step 1: T1 -> %.1f deg C  (Ys1%+.0f)    T2 -> %.1f deg C  (Ys2%+.0f)\n', ...
    Ys(1)+delta_r_step1(1), delta_r_step1(1), ...
    Ys(2)+delta_r_step1(2), delta_r_step1(2));
fprintf('  Step 2: T1 -> %.1f deg C  (Ys1%+.0f)    T2 -> %.1f deg C  (Ys2%+.0f)\n\n', ...
    Ys(1)+delta_r_step2(1), delta_r_step2(1), ...
    Ys(2)+delta_r_step2(2), delta_r_step2(2));

%% ── 4.  PI discrete parameters ──────────────────────────────────────────
q0_1 = TCL.PI.q0_1;   q1_1 = TCL.PI.q1_1;
q0_2 = TCL.PI.q0_2;   q1_2 = TCL.PI.q1_2;

u_H = [100; 100] - Us;
u_L = [  0;   0] - Us;

%% ── 5.  Setpoint filter  (beta_r = 0.95) ────────────────────────────────
beta_r    = 0.95;
phy_r     = beta_r * eye(n_op);
rk_f_now  = zeros(n_op, 1);    % running filtered setpoint (perturbation)

%% ── 6.  Innovation / error filter  (alfa_e = 0.95) ──────────────────────
alfa_e    = 0.95;
phy_e     = alfa_e * eye(n_op);
ek_f_now  = zeros(n_op, 1);    % running filtered error

%% ── 7.  PI memory ────────────────────────────────────────────────────────
e_prev = zeros(n_ip, 1);
u_prev = zeros(n_ip, 1);

%% ── 8.  Data logging ────────────────────────────────────────────────────
t1s = [];  t2s = [];
h1s = [];  h2s = [];
R1s = [];  R2s = [];
e1s = [];  e2s = [];

%% ── 9.  Live figure ──────────────────────────────────────────────────────
figure(1)

%% ── 10.  Main real-time loop ─────────────────────────────────────────────
fprintf('Starting experiment.  Total time: %.1f min\n', N_samp*Ts/60);
fprintf('Phase 1 (warm-up): heaters OFF for %d samples (%d s)\n\n', ...
    k_warmup, k_warmup*Ts);

for k = 1:N_samp
    tic

    % ---- Read hardware temperatures --------------------------------------
    Tk    = zeros(n_op,1);
    Tk(1) = T1C();
    Tk(2) = T2C();

    yk = Tk - Ys;    % perturbation from steady state

    % ---- Setpoint schedule ----------------------------------------------
    if k < k_step1
        setpt = zeros(n_op,1);         % Phases 1 & 2: hold at Ys
    elseif k < k_step2
        setpt = delta_r_step1;          % Phase 3: Step 1
    else
        setpt = delta_r_step2;          % Phase 4: Step 2
    end

    % ---- Setpoint filter ------------------------------------------------
    rk_f_now = phy_r * rk_f_now + (eye(n_op) - phy_r) * setpt;

    % ---- Error and innovation filter ------------------------------------
    raw_ek   = rk_f_now - yk;
    ek_f_now = phy_e * ek_f_now + (eye(n_op) - phy_e) * raw_ek;

    % ---- Phase 1: warm-up with heaters OFF ------------------------------
    if k <= k_warmup
        uk = [0; 0];
        % Reset PI memory so integrator starts clean
        e_prev = zeros(n_ip,1);
        u_prev = zeros(n_ip,1);
    else
        % ---- PI velocity law (Tustin) ----------------------------------
        du1 = q0_1 * ek_f_now(1) + q1_1 * e_prev(1);
        du2 = q0_2 * ek_f_now(2) + q1_2 * e_prev(2);

        u1 = u_prev(1) + du1;
        u2 = u_prev(2) + du2;

        % ---- Anti-windup clamping --------------------------------------
        u1 = min(max(u1, u_L(1)), u_H(1));
        u2 = min(max(u2, u_L(2)), u_H(2));

        uk = [u1; u2];

        % ---- Update PI memory ------------------------------------------
        e_prev = ek_f_now;
        u_prev = uk;
    end

    % ---- Absolute heater values and filtered setpoint -------------------
    Uk = uk + Us;
    Rk = rk_f_now + Ys;

    % ---- Send to heaters ------------------------------------------------
    h1(Uk(1));
    h2(Uk(2));

    % ---- Log ------------------------------------------------------------
    t1s = [t1s, Tk(1)];
    t2s = [t2s, Tk(2)];
    h1s = [h1s, Uk(1)];
    h2s = [h2s, Uk(2)];
    R1s = [R1s, Rk(1)];
    R2s = [R2s, Rk(2)];
    e1s = [e1s, ek_f_now(1)];
    e2s = [e2s, ek_f_now(2)];

    % ---- Phase banner on transitions -----------------------------------
    if k == k_warmup + 1
        fprintf('Phase 2: PI active, settling at steady state ...\n');
    elseif k == k_step1
        fprintf('Phase 3: Step 1 applied  (T1 %+.0f deg C, T2 %+.0f deg C)\n', ...
            delta_r_step1(1), delta_r_step1(2));
    elseif k == k_step2
        fprintf('Phase 4: Step 2 applied  (T1 %+.0f deg C, T2 %+.0f deg C)\n', ...
            delta_r_step2(1), delta_r_step2(2));
    end

    % ---- Live strip chart (every 15 samples = 60 s) --------------------
    if rem(k, 15) == 0
        n    = length(t1s);
        time = (0 : n-1) * Ts;

        clf
        subplot(2,1,1)
        plot(time, t1s, 'r-',  time, R1s, 'r--', ...
             time, t2s, 'b-',  time, R2s, 'b--')
        ylabel('Temperature (deg C)')
        legend('T_1','Ref_1','T_2','Ref_2','Location','NorthWest')
        title(sprintf('TCL PI Servo Experiment  –  k = %d / %d', k, N_samp))
        grid on

        subplot(2,1,2)
        plot(time, h1s, 'r-', 'LineWidth', 2)
        hold on
        plot(time, h2s, 'b-', 'LineWidth', 2)
        yline(5,  'k--', 'LineWidth', 0.8)
        yline(80, 'k--', 'LineWidth', 0.8)
        ylabel('Heater (%)')
        xlabel('Time (s)')
        legend('H_1','H_2','Location','NorthWest')
        grid on
        drawnow
    end

    % ---- Real-time pacing -----------------------------------------------
    elapsed = toc;
    fprintf('k = %3d  |  T1 = %.2f  T2 = %.2f  |  H1 = %.1f%%  H2 = %.1f%%  |  slack = %.3f s\n', ...
        k, Tk(1), Tk(2), Uk(1), Uk(2), Ts - elapsed);
    pause(max(0.01, Ts - elapsed))

end  % end main loop

%% ── 11.  Turn heaters off ────────────────────────────────────────────────
disp('Turning heaters off ...')
h1(0);
h2(0);
disp('Servo experiment complete.')

%% ── 12.  Final time axis ─────────────────────────────────────────────────
n    = length(t1s);
time = (0 : n-1) * Ts;
tf   = (N_samp - 1) * Ts;

%% ── 13.  Figure 2 – Temperature profiles ────────────────────────────────
figure(2)
subplot(2,1,1)
plot(time, t1s, 'r-', time, R1s, 'r--', 'LineWidth', 1.5), grid on
axis([0 tf min([t1s t2s])-3  max([t1s t2s])+3])
ylabel('Temperature (deg C)')
legend('T_1','Ref_1','Location','NorthWest')
title('PI Servo Experiment – Temperature Profiles')
xline(k_step1*Ts, 'k:', 'Step 1','LabelVerticalAlignment','bottom')
xline(k_step2*Ts, 'k:', 'Step 2','LabelVerticalAlignment','bottom')

subplot(2,1,2)
plot(time, t2s, 'b-', time, R2s, 'b--', 'LineWidth', 1.5), grid on
axis([0 tf min([t1s t2s])-3  max([t1s t2s])+3])
legend('T_2','Ref_2','Location','NorthWest')
ylabel('Temperature (deg C)')
xlabel('Time (s)')
xline(k_step1*Ts, 'k:')
xline(k_step2*Ts, 'k:')

%% ── 14.  Figure 3 – Heater profiles ─────────────────────────────────────
figure(3)
subplot(2,1,1)
stairs(time, h1s, 'r-', 'LineWidth', 2), grid on
axis([0 tf 0 100])
yline(5,  'k--', 'LineWidth', 0.8)
yline(80, 'k--', 'LineWidth', 0.8)
ylabel('Heater 1 (%)')
title('PI Servo Experiment – Heater Inputs')
xline(k_step1*Ts, 'k:')
xline(k_step2*Ts, 'k:')

subplot(2,1,2)
stairs(time, h2s, 'b-', 'LineWidth', 2), grid on
axis([0 tf 0 100])
yline(5,  'k--', 'LineWidth', 0.8)
yline(80, 'k--', 'LineWidth', 0.8)
ylabel('Heater 2 (%)')
xlabel('Time (s)')
xline(k_step1*Ts, 'k:')
xline(k_step2*Ts, 'k:')

%% ── 15.  Figure 4 – Filtered error ──────────────────────────────────────
figure(4)
subplot(2,1,1)
plot(time, e1s, 'r-', 'LineWidth', 1.5), grid on
axis([0 tf -8 8])
title('Filtered Error Signal  e_{f}(k)')
ylabel('e_{f,1}(k)  (deg C)')
xline(k_step1*Ts, 'k:')
xline(k_step2*Ts, 'k:')

subplot(2,1,2)
plot(time, e2s, 'b-', 'LineWidth', 1.5), grid on
axis([0 tf -8 8])
ylabel('e_{f,2}(k)  (deg C)')
xlabel('Time (s)')
xline(k_step1*Ts, 'k:')
xline(k_step2*Ts, 'k:')

%% ── 16.  Save results ────────────────────────────────────────────────────
save TCL_PI_Servo_ExpResults  time h1s h2s t1s t2s R1s R2s e1s e2s
fprintf('\nResults saved to TCL_PI_Servo_ExpResults.mat\n');
