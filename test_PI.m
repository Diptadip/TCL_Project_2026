% =========================================================================
%  TCL_PI_Experiment.m
%
%  Runs a real-time experiment on the TCL Arduino board using two
%  decentralised IMC-tuned PI controllers (H1->T1, H2->T2).
%
%  PI controller – velocity (incremental) form, Tustin discretisation:
%      u(k) = u(k-1) + q0*e(k) + q1*e(k-1)
%
%  Setpoint schedule (mirrors LQOC experimental script):
%    k =   1 – 49   : manual mode  (heaters OFF, system at rest)
%    k =  50 – 249  : PI active, setpoint = steady state  (no step yet)
%    k = 250 – 449  : PI active, setpoint step  +8 deg C on T1, +5 on T2
%
%  Plots produced (live + post-experiment, matching LQOC script):
%    Fig 1 – Live temperature and heater strip chart (updated every 15 s)
%    Fig 2 – Temperature profiles vs setpoint (post-experiment)
%    Fig 3 – Heater input profiles            (post-experiment)
%    Fig 4 – Filtered error / mismatch signal (post-experiment)
%
%  Dependencies:
%    tclab_N.m                     (Arduino initialisation)
%    TCL_MechModel_Parameters.mat  (contains TCL struct with PI sub-struct)
% =========================================================================

close all; clear all; clc

%% ── 0.  Connect to Arduino TCL board ────────────────────────────────────
tclab_N;

%% ── 1.  Load PI parameters ───────────────────────────────────────────────
N_samp = 450;
disp('Load PI Controller Parameters')

load TCL_MechModel_Parameters.mat   % loads struct TCL

if ~isfield(TCL, 'PI')
    error(['TCL.PI not found. Run TCL_PI_Controller_Design.m first ', ...
           'to compute and save PI tuning parameters.']);
end

Ys = TCL.Ys;        % steady-state output  (column vector, deg C)
Us = TCL.Us;        % steady-state input   (column vector, %)
Ts = TCL.Samp_T;    % sampling period (s)

n_ip = 2;
n_op = 2;

%% ── 2.  PI discrete parameters ──────────────────────────────────────────
q0_1 = TCL.PI.q0_1;    q1_1 = TCL.PI.q1_1;   % Loop 1: H1 -> T1
q0_2 = TCL.PI.q0_2;    q1_2 = TCL.PI.q1_2;   % Loop 2: H2 -> T2

%  Physical heater bounds (perturbation form)
u_H = [100; 100] - Us;
u_L = [  0;   0] - Us;

%% ── 3.  Setpoint filter  (beta_r = 0.95, same as simulation) ────────────
beta_r = 0.95;
phy_r  = beta_r * eye(n_op);

% Running filtered setpoint (column vector, perturbation form)
rk_f_now = zeros(n_op, 1);

%% ── 4.  Innovation / error filter  (alfa_e = 0.95) ──────────────────────
alfa_e   = 0.95;
phy_e    = alfa_e * eye(n_op);

ek_f_now = zeros(n_op, 1);   % running filtered error

%% ── 5.  PI controller memory ─────────────────────────────────────────────
e_prev = zeros(n_ip, 1);   % e(k-1)
u_prev = zeros(n_ip, 1);   % u(k-1)  (perturbation)

%% ── 6.  Data logging arrays (same names as LQOC script) ─────────────────
t1s = [];  t2s = [];    % measured temperatures
h1s = [];  h2s = [];    % absolute heater outputs (%)
R1s = [];  R2s = [];    % absolute filtered setpoints (deg C)
e1s = [];  e2s = [];    % filtered errors

%% ── 7.  Open live plot figure ────────────────────────────────────────────
figure(1)

%% ── 8.  Main real-time loop ──────────────────────────────────────────────
for k = 1:N_samp
    tic

    % ---- Read temperatures from TCL board --------------------------------
    Tk    = zeros(n_op, 1);
    Tk(1) = T1C();
    Tk(2) = T2C();

    yk = Tk - Ys;           % output perturbation from steady state

    % ---- Setpoint schedule -----------------------------------------------
    if k < 250
        setpt = zeros(n_op, 1);     % hold at steady state
    else
        setpt = [8; 5];             % +8 deg C on T1, +5 deg C on T2
    end

    % ---- Setpoint filter -------------------------------------------------
    rk_f_now = phy_r * rk_f_now + (eye(n_op) - phy_r) * setpt;

    % ---- Raw error -------------------------------------------------------
    raw_ek = rk_f_now - yk;

    % ---- Innovation (error) filter ---------------------------------------
    ek_f_now = phy_e * ek_f_now + (eye(n_op) - phy_e) * raw_ek;

    % ---- Manual mode: heaters OFF for first 50 samples ------------------
    if k <= 50
        uk = [0; 0];
    else
        % ---- PI velocity law -------------------------------------------
        %  Loop 1  (H1 -> T1)
        du1 = q0_1 * ek_f_now(1) + q1_1 * e_prev(1);
        u1  = u_prev(1) + du1;

        %  Loop 2  (H2 -> T2)
        du2 = q0_2 * ek_f_now(2) + q1_2 * e_prev(2);
        u2  = u_prev(2) + du2;

        % ---- Anti-windup clamping --------------------------------------
        u1 = min(max(u1, u_L(1)), u_H(1));
        u2 = min(max(u2, u_L(2)), u_H(2));

        uk = [u1; u2];

        % ---- Update PI memory ------------------------------------------
        e_prev = ek_f_now;
        u_prev = uk;
    end

    % ---- Absolute heater and setpoint values ----------------------------
    Uk = uk + Us;                        % absolute heater power (%)
    Rk = rk_f_now + Ys;                  % absolute filtered setpoint (deg C)

    % ---- Send to heaters ------------------------------------------------
    h1(Uk(1));
    h2(Uk(2));

    % ---- Log data (same variable names as LQOC script) ------------------
    t1s = [t1s, Tk(1)];
    t2s = [t2s, Tk(2)];
    h1s = [h1s, Uk(1)];
    h2s = [h2s, Uk(2)];
    R1s = [R1s, Rk(1)];
    R2s = [R2s, Rk(2)];
    e1s = [e1s, ek_f_now(1)];
    e2s = [e2s, ek_f_now(2)];

    % ---- Live strip chart (update every 15 samples) ---------------------
    if rem(k, 15) == 0
        n    = length(t1s);
        time = linspace(0, n * Ts, n);

        clf
        subplot(2,1,1)
        plot(time, t1s, 'r-',  time, R1s, 'k.-', 'MarkerSize', 10)
        hold on
        plot(time, t2s, 'b-',  time, R2s, 'k.-', 'MarkerSize', 10)
        ylabel('Temperature (degC)')
        legend('Temperature 1','Temperature 2','Location','NorthWest')
        title('TCL – PI Control (live)')

        subplot(2,1,2)
        plot(time, h1s, 'r-', 'LineWidth', 2)
        hold on
        plot(time, h2s, 'b-', 'LineWidth', 2)
        ylabel('Heater (%)')
        xlabel('Time (sec)')
        legend('Heater 1','Heater 2','Location','NorthWest')
        drawnow
    end

    % ---- Real-time pacing -----------------------------------------------
    elapsed = toc;
    fprintf('k = %d   time remaining in sample = %.3f s\n', k, Ts - elapsed)
    pause(max(0.01, Ts - elapsed))

end  % end main loop

%% ── 9.  Turn off heaters ─────────────────────────────────────────────────
disp('Turning off heaters')
h1(0);
h2(0);
disp('PI Experimental Test Complete')

%% ── 10.  Build final time axis ───────────────────────────────────────────
tf   = N_samp * Ts;
n    = length(t1s);
time = linspace(0, n * Ts, n);

%% ── 11.  Figure 2 – Temperature profiles (matches LQOC Fig 2) ───────────
figure(2)
subplot(2,1,1)
plot(time, t1s, 'r-', time, R1s, 'k.-', 'MarkerSize', 10), grid
axis([0 tf 30 70])
ylabel('Temperature (degC)')
legend('Temperature 1','Setpoint 1')
title('PI Control – Temperature Profiles')

subplot(2,1,2)
plot(time, t2s, 'b-', time, R2s, 'k.-', 'MarkerSize', 10), grid
axis([0 tf 30 70])
legend('Temperature 2','Setpoint 2')
ylabel('Temperature (degC)')
xlabel('Time (sec)')

%% ── 12.  Figure 3 – Heater input profiles (matches LQOC Fig 3) ──────────
figure(3)
subplot(2,1,1)
stairs(time, h1s, 'r-', 'LineWidth', 2), grid
axis([0 tf 0 100])
ylabel('Heater 1 (%)')
title('PI Control – Heater Input Profiles')

subplot(2,1,2)
stairs(time, h2s, 'b-', 'LineWidth', 2), grid
axis([0 tf 0 100])
ylabel('Heater 2 (%)')
xlabel('Time (sec)')

%% ── 13.  Figure 4 – Filtered error signal (matches LQOC Fig 4) ──────────
figure(4)
subplot(2,1,1)
plot(time, e1s, 'r-', 'MarkerSize', 10), grid
axis([0 tf -5 5])
title('Filtered Error Signal  e_{f}(k)')
ylabel('e_{f,1}(k)  (degC)')

subplot(2,1,2)
plot(time, e2s, 'b-', 'MarkerSize', 10), grid
axis([0 tf -5 5])
ylabel('e_{f,2}(k)  (degC)')
xlabel('Time (sec)')

%% ── 14.  Save experimental results ──────────────────────────────────────
save TCL_PI_Experiment_Results  time h1s h2s t1s t2s R1s R2s e1s e2s
disp('Results saved to TCL_PI_Experiment_Results.mat')