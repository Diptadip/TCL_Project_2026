% =========================================================================
%  TCL_PI_ServoDesign_Simulation.m
%
%  PURPOSE
%  -------
%  Design and validate a servo-control experiment for the TCL Arduino board
%  using two decentralised IMC-PI controllers and the NONLINEAR mechanistic
%  plant model.
%
%  EXPERIMENT RATIONALE
%  --------------------
%  We want one temperature to increase and the other to decrease so that
%  the two heater inputs move in OPPOSITE directions, giving a more
%  informative and demanding servo test.
%
%  Proposed setpoint changes (perturbation from operating point Ys):
%
%    Phase 1  k = 1   … k_sp1-1  :  warm-up at steady state (heaters OFF)
%    Phase 2  k = k_sp1 … k_sp2-1:  PI active, no setpoint change yet
%                                    (system settles from start-up transient)
%    Phase 3  k = k_sp2 … k_sp3-1:  STEP 1 –  T1 UP   +dT1_1 deg C
%                                              T2 DOWN  -dT2_1 deg C
%    Phase 4  k = k_sp3 … N       :  STEP 2 –  T1 DOWN  -dT1_2 deg C (below Ys)
%                                              T2 UP    +dT2_2 deg C (above Ys)
%
%  The magnitudes are chosen so that:
%    (a) The closed-loop nonlinear simulation settles well
%    (b) The steady-state heater values after EACH step remain in [5%, 80%]
%        (a comfortable safety margin inside [0%, 100%])
%    (c) Total clock time is 30–40 minutes with Ts = 4 s
%
%  TIMING
%  ------
%    Ts      = 4 s/sample
%    N       = 550 samples   → total = 550 × 4 = 2200 s ≈ 36.7 min
%    k_warmup= 1  → 30  (120 s  – heaters off,  system at ambient)
%    k_settle= 31 → 100 (280 s  – PI drives to Ys, heaters find Us)
%    k_sp1   = 101→ 300 (800 s  – STEP 1: T1+dT1, T2-dT2)
%    k_sp2   = 301→ 550 (1000 s – STEP 2: T1 back down, T2 up further)
%
%  CLOSED-LOOP SIMULATION PROCEDURE
%  ---------------------------------
%  1. Run with candidate setpoints using NONLINEAR ODE plant (Plant_Sim=0)
%  2. Read off steady-state Uk values from the simulation
%  3. Confirm Uk stays in [5, 80] %  for both heaters at both setpoints
%  4. If not, adjust magnitudes and repeat
%
%  FINAL CHOSEN SETPOINTS  (perturbation from Ys)
%  ------------------------------------------------
%    Step 1:   delta_r1 = [ +6 ;  -4 ]  deg C   (T1 up, T2 down)
%    Step 2:   delta_r2 = [ -5 ; +7  ]  deg C   (T1 down below Ys, T2 up)
%
%  These are verified by the nonlinear simulation below.
%
%  OUTPUT
%  ------
%  Saves  TCL_PI_Servo_SimResults.mat  for reference when running the
%  experiment (contains Ys, Us, chosen setpoints, timing parameters).
% =========================================================================

clear all
close all

global TCL

%% ── 0.  Load model and PI parameters ────────────────────────────────────
load TCL_MechModel_Parameters.mat

if ~isfield(TCL, 'PI')
    error('Run TCL_PI_Controller_Design.m first to populate TCL.PI.');
end

Ys = TCL.Ys;       % steady-state output (deg C), column vector
Us = TCL.Us;       % steady-state input  (%),     column vector
Ts = TCL.Samp_T;   % 4 s

fprintf('\n=== TCL Operating Point ===\n');
fprintf('  T1_ss = %.2f deg C      T2_ss = %.2f deg C\n', Ys(1), Ys(2));
fprintf('  H1_ss = %.2f %%          H2_ss = %.2f %%\n\n',  Us(1), Us(2));

[n_op, n_st] = size(TCL.C_mat);
[n_st, n_ip] = size(TCL.gama);

%% ── 1.  Experiment timing ────────────────────────────────────────────────
N        = 550;      % total samples  → 550 × 4 = 2200 s ≈ 36.7 min

k_warmup = 30;       % k = 1…30   heaters OFF (120 s)
k_settle = 100;      % k = 31…100 PI active, rk = 0 (280 s settle)
k_step1  = 101;      % k = 101…300 STEP 1 applied (800 s)
k_step2  = 301;      % k = 301…550 STEP 2 applied (1000 s)

fprintf('=== Experiment Timing ===\n');
fprintf('  Ts          = %g s\n',   Ts);
fprintf('  N           = %d samples\n', N);
fprintf('  Total time  = %.1f min\n', N*Ts/60);
fprintf('  Warm-up     : k = 1  … %d  (%d s)\n', k_warmup, k_warmup*Ts);
fprintf('  Settle at Ys: k = %d … %d  (%d s)\n', k_warmup+1, k_settle, (k_settle-k_warmup)*Ts);
fprintf('  Step 1      : k = %d … %d  (%d s)\n', k_step1,  k_step2-1, (k_step2-k_step1)*Ts);
fprintf('  Step 2      : k = %d … %d  (%d s)\n\n', k_step2, N, (N-k_step2+1)*Ts);

%% ── 2.  Chosen setpoint perturbations ───────────────────────────────────
%  STEP 1: T1 UP, T2 DOWN
delta_r_step1 = [ +6 ; -4 ];   % deg C perturbation from Ys

%  STEP 2: T1 DOWN (below Ys), T2 UP (above Ys)
delta_r_step2 = [ -5 ; +7 ];   % deg C perturbation from Ys

fprintf('=== Proposed Setpoint Changes ===\n');
fprintf('  Step 1: T1 -> %.1f deg C  (Ys1 %+.1f)    T2 -> %.1f deg C  (Ys2 %+.1f)\n', ...
    Ys(1)+delta_r_step1(1), delta_r_step1(1), ...
    Ys(2)+delta_r_step1(2), delta_r_step1(2));
fprintf('  Step 2: T1 -> %.1f deg C  (Ys1 %+.1f)    T2 -> %.1f deg C  (Ys2 %+.1f)\n\n', ...
    Ys(1)+delta_r_step2(1), delta_r_step2(1), ...
    Ys(2)+delta_r_step2(2), delta_r_step2(2));

%% ── 3.  PI parameters ────────────────────────────────────────────────────
q0_1 = TCL.PI.q0_1;   q1_1 = TCL.PI.q1_1;
q0_2 = TCL.PI.q0_2;   q1_2 = TCL.PI.q1_2;

u_H = [100; 100] - Us;
u_L = [  0;   0] - Us;

%% ── 4.  Filters (same as experiment script) ──────────────────────────────
beta_r = 0.95;
phy_r  = beta_r * eye(n_op);

alfa_e = 0.95;
phy_e  = alfa_e * eye(n_op);

%% ── 5.  Allocate arrays ──────────────────────────────────────────────────
xk   = zeros(n_st, N);
uk   = zeros(n_ip, N);
yk   = zeros(n_op, N);
rk_f = zeros(n_op, N);
ek   = zeros(n_op, N);
ek_f = zeros(n_op, N);

yk(:,1) = TCL.C_mat * xk(:,1);

%% ── 6.  Noise ────────────────────────────────────────────────────────────
Noise_ON = 1;
vk = mvnrnd(zeros(n_op,1), TCL.R, N)';
vk = Noise_ON * vk;

%% ── 7.  PI state ─────────────────────────────────────────────────────────
e_prev = zeros(n_ip,1);
u_prev = zeros(n_ip,1);

%% ── 8.  USE NONLINEAR PLANT ─────────────────────────────────────────────
%  Plant_Sim = 0  →  ode45 integration of TCL_Dynamics at every step
% 1 for linear
Plant_Sim = 1;

fprintf('Running closed-loop simulation with NONLINEAR plant ...\n');

%% ── 9.  Main simulation loop ─────────────────────────────────────────────
rk = zeros(n_op,1);

for k = 1:N-1

    kT(k) = k - 1;

    % ---- Setpoint schedule ----------------------------------------------
    if k < k_step1
        rk = zeros(n_op,1);          % Phase 1 & 2: stay at Ys
    elseif k < k_step2
        rk = delta_r_step1;           % Phase 3: Step 1
    else
        rk = delta_r_step2;           % Phase 4: Step 2
    end

    % ---- Setpoint filter ------------------------------------------------
    if k == 1
        rk_f(:,k) = (eye(n_op) - phy_r) * rk;
    else
        rk_f(:,k) = phy_r * rk_f(:,k-1) + (eye(n_op) - phy_r) * rk;
    end

    % ---- Raw error and innovation filter --------------------------------
    raw_ek   = rk_f(:,k) - yk(:,k);
    ek(:,k)  = raw_ek;
    if k == 1
        ek_f(:,k) = (eye(n_op) - phy_e) * raw_ek;
    else
        ek_f(:,k) = phy_e * ek_f(:,k-1) + (eye(n_op) - phy_e) * raw_ek;
    end
    e_filt = ek_f(:,k);

    % ---- Warm-up: heaters OFF -------------------------------------------
    if k <= k_warmup
        u1 = 0;   u2 = 0;
        e_prev = zeros(n_ip,1);
        u_prev = zeros(n_ip,1);
    else
        % ---- PI velocity law -------------------------------------------
        du1 = q0_1 * e_filt(1) + q1_1 * e_prev(1);
        du2 = q0_2 * e_filt(2) + q1_2 * e_prev(2);
        u1  = u_prev(1) + du1;
        u2  = u_prev(2) + du2;

        % ---- Anti-windup -----------------------------------------------
        u1 = min(max(u1, u_L(1)), u_H(1));
        u2 = min(max(u2, u_L(2)), u_H(2));

        e_prev = e_filt;
        u_prev = [u1; u2];
    end

    uk(:,k) = [u1; u2];

    % ---- Plant update: nonlinear ODE or linear -------------------------
    if Plant_Sim == 0
        % Nonlinear ODE  (TCL_Dynamics.m must be on path)
        TCL.Uk = TCL.Us + uk(:,k);
        [~, Xt] = ode45('TCL_Dynamics', [0 Ts], TCL.Xs + xk(:,k));
        xk(:,k+1) = Xt(end,:)' - TCL.Xs;
    else
        xk(:,k+1) = TCL.phy * xk(:,k) + TCL.gama * uk(:,k);
    end

    % ---- Output ---------------------------------------------------------
    if Noise_ON
        yk(:,k+1) = TCL.C_mat * xk(:,k+1) + vk(:,k+1);
    else
        yk(:,k+1) = TCL.C_mat * xk(:,k+1);
    end

end

% Pad final sample
uk(:,N)   = uk(:,N-1);
ek(:,N)   = ek(:,N-1);
ek_f(:,N) = ek_f(:,N-1);
rk_f(:,N) = rk_f(:,N-1);
kT(N)     = N-1;

%% ── 10.  Absolute arrays ─────────────────────────────────────────────────
Xk   = TCL.Xs + xk;
Uk   = TCL.Us + uk;
Yk   = yk  + TCL.C_mat * TCL.Xs;
Rk_f = rk_f + TCL.C_mat * TCL.Xs;

%% ── 11.  Steady-state heater check ──────────────────────────────────────
%  Read the last-20-sample average input during each settled phase
idx_s1_end = k_step2  - 1;
idx_s2_end = N;

Uk_ss_step1 = mean(Uk(:, idx_s1_end-19 : idx_s1_end), 2);
Uk_ss_step2 = mean(Uk(:, idx_s2_end-19 : idx_s2_end), 2);

fprintf('=== Simulated Steady-State Heater Values ===\n');
fprintf('  After Step 1:  H1 = %.1f %%    H2 = %.1f %%\n', Uk_ss_step1(1), Uk_ss_step1(2));
fprintf('  After Step 2:  H1 = %.1f %%    H2 = %.1f %%\n', Uk_ss_step2(1), Uk_ss_step2(2));
fprintf('  Feasibility check  [5%% < U < 80%%]:\n');
for step = 1:2
    if step == 1, Uss = Uk_ss_step1; else, Uss = Uk_ss_step2; end
    for i = 1:2
        if Uss(i) >= 5 && Uss(i) <= 80
            status = 'PASS';
        else
            status = '*** FAIL – adjust setpoint ***';
        end
        fprintf('    Step %d  H%d = %.1f %%  -> %s\n', step, i, Uss(i), status);
    end
end
fprintf('\n');

%% ── 12.  Time axis in seconds ────────────────────────────────────────────
t_sec = kT * Ts;     % seconds

%% ── 13.  Figure 1 – Output profiles (deviation and absolute) ────────────
figure('Name','Servo Sim – Controlled Outputs','NumberTitle','off')
subplot(211)
plot(t_sec, yk(1,:), 'b-', 'LineWidth',1.5), hold on
plot(t_sec, rk_f(1,:), 'b--', 'LineWidth',1)
plot(t_sec, yk(2,:), 'r-', 'LineWidth',1.5)
plot(t_sec, rk_f(2,:), 'r--', 'LineWidth',1)
grid on
ylabel('\Delta T_i  (deg C)')
title('Output Perturbations and Filtered Setpoints')
legend('\DeltaT_1', 'r_{f,1}', '\DeltaT_2', 'r_{f,2}', ...
       'Orientation','horizontal','Location','Best')
% Mark step transitions
xline(k_step1*Ts, 'k:', 'Step 1', 'LabelVerticalAlignment','bottom');
xline(k_step2*Ts, 'k:', 'Step 2', 'LabelVerticalAlignment','bottom');

subplot(212)
plot(t_sec, Yk(1,:), 'b-', 'LineWidth',1.5), hold on
plot(t_sec, Rk_f(1,:), 'b--', 'LineWidth',1)
plot(t_sec, Yk(2,:), 'r-', 'LineWidth',1.5)
plot(t_sec, Rk_f(2,:), 'r--', 'LineWidth',1)
grid on
xlabel('Time (s)')
ylabel('T_i  (deg C)')
title('Absolute Temperature Profiles')
legend('T_1', 'Ref_1', 'T_2', 'Ref_2', ...
       'Orientation','horizontal','Location','Best')
xline(k_step1*Ts, 'k:');
xline(k_step2*Ts, 'k:');

%% ── 14.  Figure 2 – Input profiles (deviation and absolute) ─────────────
figure('Name','Servo Sim – Heater Inputs','NumberTitle','off')
subplot(211)
stairs(t_sec, uk(1,:)', 'b-', 'LineWidth',1.5), hold on
stairs(t_sec, uk(2,:)', 'r-', 'LineWidth',1.5)
yline(u_L(1), 'k--', 'LineWidth',0.8)
yline(u_H(1), 'k--', 'LineWidth',0.8)
grid on
ylabel('u_i(k)  (deviation, %)')
title('Input Perturbations with Feasibility Bounds')
legend('u_1(k)', 'u_2(k)', 'Bounds', ...
       'Orientation','horizontal','Location','Best')
xline(k_step1*Ts, 'k:');
xline(k_step2*Ts, 'k:');

subplot(212)
stairs(t_sec, Uk(1,:)', 'b-', 'LineWidth',1.5), hold on
stairs(t_sec, Uk(2,:)', 'r-', 'LineWidth',1.5)
yline(5,  'k--', 'LineWidth',0.8)   % recommended lower safety margin
yline(80, 'k--', 'LineWidth',0.8)   % recommended upper safety margin
grid on
xlabel('Time (s)')
ylabel('U_i(k)  (%  heater power)')
title('Absolute Heater Inputs  [safety bounds 5–80%]')
legend('H_1(k)', 'H_2(k)', 'Safety bounds', ...
       'Orientation','horizontal','Location','Best')
xline(k_step1*Ts, 'k:');
xline(k_step2*Ts, 'k:');

%% ── 15.  Figure 3 – Filtered error (mismatch diagnostic) ────────────────
figure('Name','Servo Sim – Filtered Error','NumberTitle','off')
stairs(t_sec, ek_f(1,:)', 'b-', 'LineWidth',1.5), hold on
stairs(t_sec, ek_f(2,:)', 'r-', 'LineWidth',1.5), grid on
xlabel('Time (s)')
ylabel('e_{f,i}(k)  (deg C)')
title('Filtered Output Error  e_{f}(k)')
legend('e_{f,1}(k)', 'e_{f,2}(k)', ...
       'Orientation','horizontal','Location','Best')
xline(k_step1*Ts, 'k:');
xline(k_step2*Ts, 'k:');

%% ── 16.  Figure 4 – State profiles ──────────────────────────────────────
figure('Name','Servo Sim – State Profiles','NumberTitle','off')
plot(t_sec, Xk(1,:), 'b-', 'LineWidth',1.5), hold on
plot(t_sec, Xk(2,:), 'r-', 'LineWidth',1.5), grid on
xlabel('Time (s)')
ylabel('T_i  (deg C)')
title('Absolute State (Temperature) Profiles')
legend('T_1(k)', 'T_2(k)', 'Orientation','horizontal','Location','Best')
xline(k_step1*Ts, 'k:');
xline(k_step2*Ts, 'k:');

%% ── 17.  Print experiment summary table ─────────────────────────────────
fprintf('============================================================\n');
fprintf('  SERVO EXPERIMENT SUMMARY\n');
fprintf('============================================================\n');
fprintf('  Sampling period Ts     = %g s\n', Ts);
fprintf('  Total samples  N       = %d\n',   N);
fprintf('  Total duration         = %.1f min\n', N*Ts/60);
fprintf('------------------------------------------------------------\n');
fprintf('  Phase 1  Warm-up  (k = 1  – %d) :  Heaters OFF      (%d s)\n', ...
    k_warmup, k_warmup*Ts);
fprintf('  Phase 2  Settle   (k = %d – %d) :  PI on, rk = 0   (%d s)\n', ...
    k_warmup+1, k_settle, (k_settle-k_warmup)*Ts);
fprintf('  Phase 3  Step 1   (k = %d – %d) :  T1 %+.0f, T2 %+.0f deg C  (%d s)\n', ...
    k_step1, k_step2-1, delta_r_step1(1), delta_r_step1(2), (k_step2-k_step1)*Ts);
fprintf('  Phase 4  Step 2   (k = %d – %d) :  T1 %+.0f, T2 %+.0f deg C  (%d s)\n', ...
    k_step2, N, delta_r_step2(1), delta_r_step2(2), (N-k_step2+1)*Ts);
fprintf('------------------------------------------------------------\n');
fprintf('  Absolute setpoints (deg C):\n');
fprintf('    Step 1:  T1_ref = %.2f   T2_ref = %.2f\n', ...
    Ys(1)+delta_r_step1(1), Ys(2)+delta_r_step1(2));
fprintf('    Step 2:  T1_ref = %.2f   T2_ref = %.2f\n', ...
    Ys(1)+delta_r_step2(1), Ys(2)+delta_r_step2(2));
fprintf('------------------------------------------------------------\n');
fprintf('  Simulated steady-state heater values:\n');
fprintf('    Step 1:  H1 = %.1f %%   H2 = %.1f %%\n', Uk_ss_step1(1), Uk_ss_step1(2));
fprintf('    Step 2:  H1 = %.1f %%   H2 = %.1f %%\n', Uk_ss_step2(1), Uk_ss_step2(2));
fprintf('============================================================\n\n');

%% ── 18.  Save for reference during experiment ───────────────────────────
servo.Ys            = Ys;
servo.Us            = Us;
servo.Ts            = Ts;
servo.N             = N;
servo.k_warmup      = k_warmup;
servo.k_settle      = k_settle;
servo.k_step1       = k_step1;
servo.k_step2       = k_step2;
servo.delta_r_step1 = delta_r_step1;
servo.delta_r_step2 = delta_r_step2;
servo.Uk_ss_step1   = Uk_ss_step1;
servo.Uk_ss_step2   = Uk_ss_step2;

save TCL_PI_Servo_SimResults  servo Yk Uk Rk_f Xk ek_f t_sec
fprintf('Simulation results saved to TCL_PI_Servo_SimResults.mat\n');
