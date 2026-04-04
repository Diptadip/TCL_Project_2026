% =========================================================================
%  TCL_PI_Simulation.m
%
%  Simulates the TCL closed-loop system under two decentralised
%  IMC-tuned PI controllers and produces the same set of plots as the
%  LQOC simulation script (TCL_LQOC_MechMod_Simulation.m).
%
%  PI controller – velocity (incremental) form, Tustin discretisation:
%      u(k) = u(k-1) + q0*e(k) + q1*e(k-1)
%
%  Plots produced (matching LQOC script):
%    Fig 1 – Controlled output and setpoint profiles (deviation + absolute)
%    Fig 2 – Manipulated input profiles (deviation + absolute)
%    Fig 3 – Model-plant mismatch signal  e_f(k)  [innovation / error]
%    Fig 4 – System state profiles
%    Fig 5 – PI internal signals: filtered setpoint and error per loop
%
%  Dependencies:
%    TCL_MechModel_Parameters.mat  (contains TCL struct with phy, gama,
%                                   C_mat, Xs, Us, Ys, Samp_T, R, PI)
%    Set_Graphics_Parameters       (graphics defaults script, same as LQOC)
% =========================================================================

clear all
close all

global TCL

%% ── 0.  Load model and PI parameters ────────────────────────────────────
load TCL_MechModel_Parameters.mat
% TCL_PI_Controller_Design.m must have been run first so that TCL.PI exists.
% If TCL.PI is missing, run that script now to populate it.
if ~isfield(TCL, 'PI')
    error(['TCL.PI not found in TCL_MechModel_Parameters.mat. ', ...
           'Run TCL_PI_Controller_Design.m first.']);
end

%% ── 1.  System dimensions and steady-state point ────────────────────────
[n_op, ~] = size(TCL.C_mat);          % n_op = 2 outputs, n_st = 2 states
[n_st, n_ip] = size(TCL.gama);           % n_ip = 2 inputs

Ys = TCL.Ys;    % steady-state output  (column vector, deg C)
Us = TCL.Us;    % steady-state input   (column vector, %)
Ts = TCL.Samp_T;

%% ── 2.  Simulation length ────────────────────────────────────────────────
N = 400;        % number of sampling instants  (same as LQOC script)

%% ── 3.  PI discrete parameters ──────────────────────────────────────────
%  Velocity form:  u(k) = u(k-1) + q0*e(k) + q1*e(k-1)
q0_1 = TCL.PI.q0_1;    q1_1 = TCL.PI.q1_1;   % Loop 1: H1 -> T1
q0_2 = TCL.PI.q0_2;    q1_2 = TCL.PI.q1_2;   % Loop 2: H2 -> T2

%  Input perturbation bounds (same anti-windup limits as LQOC script)
u_H = [100 100]' - Us;   % upper perturbation limit
u_L = [0   0  ]' - Us;   % lower perturbation limit

%% ── 4.  Setpoint filter  (mirrors LQR.phy_r in LQOC script) ─────────────
%  A first-order exponential filter on the reference smooths the step
%  change so the heater input does not jump immediately.
%      rk_f(k) = beta_r * rk_f(k-1) + (1-beta_r) * rk(k)
beta_r  = 0.95;                          % same default as LQOC script
phy_r   = beta_r * eye(n_op);
rk_f    = zeros(n_op, N);               % filtered setpoint history

%% ── 5.  Innovation / error filter  (mirrors LQR.phy_e) ──────────────────
%  Low-pass filters the raw output error to suppress measurement noise
%  before it enters the PI difference equations.
%      ek_f(k) = alfa_e * ek_f(k-1) + (1-alfa_e) * ek(k)
alfa_e  = 0.95;                          % same default as LQOC script
phy_e   = alfa_e * eye(n_op);
ek      = zeros(n_op, N);               % raw error history
ek_f    = zeros(n_op, N);              % filtered error history

%% ── 6.  Pre-allocate state, input, output arrays ────────────────────────
xk = zeros(n_st, N);    % state perturbation
uk = zeros(n_ip, N);    % input perturbation
yk = zeros(n_op, N);    % output perturbation (= measurement - Ys)

yk(:,1) = TCL.C_mat * xk(:,1);   % initial output perturbation = 0

%% ── 7.  Measurement noise ────────────────────────────────────────────────
Noise_ON = 1;
vk = mvnrnd(zeros(n_op,1), TCL.R, N);
vk = Noise_ON * vk';      % (n_op x N)

%% ── 8.  PI controller internal state ────────────────────────────────────
e_prev = zeros(n_ip, 1);   % e(k-1) for each loop
u_prev = zeros(n_ip, 1);   % u(k-1) for each loop

%% ── 9.  Reference / setpoint (perturbation form) ────────────────────────
rk = zeros(n_op, 1);       % starts at zero (at steady state)

%% ── 10. Plant simulation choice ──────────────────────────────────────────
%  Plant_Sim = 1 : linear discrete model  (same as LQOC default)
%  Plant_Sim = 0 : nonlinear ODE  (requires TCL_Dynamics.m)
Plant_Sim = 1;

%% ── 11. Main simulation loop ─────────────────────────────────────────────
for k = 1:N-1

    kT(k) = k - 1;   % time axis in sampling instants

    % ------ Setpoint change (same schedule as LQOC script) ---------------
    if k >= 50
        rk = [8; 5];    % +8 deg C on T1, +5 deg C on T2  (perturbation)
    end

    % ------ Setpoint filter  rk_f(k) = phy_r * rk_f(k-1) + (I-phy_r)*rk
    if k == 1
        rk_f(:,k) = (eye(n_op) - phy_r) * rk;
    else
        rk_f(:,k) = phy_r * rk_f(:,k-1) + (eye(n_op) - phy_r) * rk;
    end

    % ------ Raw output error (using filtered setpoint) -------------------
    raw_ek = rk_f(:,k) - yk(:,k);
    ek(:,k) = raw_ek;

    % ------ Innovation filter  ek_f(k) = phy_e*ek_f(k-1) + (I-phy_e)*ek
    if k == 1
        ek_f(:,k) = (eye(n_op) - phy_e) * raw_ek;
    else
        ek_f(:,k) = phy_e * ek_f(:,k-1) + (eye(n_op) - phy_e) * raw_ek;
    end

    e_filt = ek_f(:,k);   % filtered error used by PI

    % ------ PI velocity law (two independent loops) ----------------------
    %  Loop 1  (H1 -> T1)
    du1 = q0_1 * e_filt(1) + q1_1 * e_prev(1);
    u1  = u_prev(1) + du1;

    %  Loop 2  (H2 -> T2)
    du2 = q0_2 * e_filt(2) + q1_2 * e_prev(2);
    u2  = u_prev(2) + du2;

    % ------ Anti-windup: clamp to physical heater range ------------------
    u1 = min(max(u1, u_L(1)), u_H(1));
    u2 = min(max(u2, u_L(2)), u_H(2));

    uk(:,k) = [u1; u2];

    % ------ Update PI memory ---------------------------------------------
    e_prev = e_filt;
    u_prev = [u1; u2];

    % ------ Plant simulation: k -> k+1 -----------------------------------
    if Plant_Sim
        % Linear discrete perturbation model
        xk(:,k+1) = TCL.phy * xk(:,k) + TCL.gama * uk(:,k);
    else
        % Nonlinear ODE simulation (requires TCL_Dynamics.m)
        TCL.Uk = TCL.Us + uk(:,k);
        [~, Xt] = ode45('TCL_Dynamics', [0 Ts], TCL.Xs + xk(:,k));
        xk(:,k+1) = Xt(end,:)' - TCL.Xs;
    end

    % ------ Output with measurement noise --------------------------------
    if Noise_ON
        yk(:,k+1) = TCL.C_mat * xk(:,k+1) + vk(:,k+1);
    else
        yk(:,k+1) = TCL.C_mat * xk(:,k+1);
    end

end  % end main loop

% ── Pad final sample (mirrors LQOC end-padding) ──────────────────────────
uk(:,N)    = uk(:,N-1);
ek(:,N)    = ek(:,N-1);
ek_f(:,N)  = ek_f(:,N-1);
rk_f(:,N)  = rk_f(:,N-1);
kT(N)      = N - 1;

%% ── 12.  Absolute (physical) variable arrays ─────────────────────────────
%  Mirrors the LQOC script's  Xk, Uk, Yk, Rk_f  calculations
Xk   = TCL.Xs + xk;                           % absolute states  (deg C)
Uk   = TCL.Us + uk;                            % absolute inputs  (%)
Yk   = yk   + TCL.C_mat * TCL.Xs;             % absolute outputs (deg C)
Rk_f = rk_f + TCL.C_mat * TCL.Xs;            % absolute filtered setpoint

u_limit = [u_L(1)*ones(N,1)  u_H(1)*ones(N,1)];   % perturbation bounds table

%% ── 13.  Graphics defaults ───────────────────────────────────────────────
Set_Graphics_Parameters    % same call as LQOC script

%% ── 14.  Figure 1 – Controlled output and setpoint profiles ─────────────
%  Structure mirrors the LQOC Figure 1 exactly.
figure
subplot(211)
plot(kT, yk(1,:),    'b',  ...
     kT, rk_f(1,:), 'b-.', ...
     kT, yk(2,:),    'k',  ...
     kT, rk_f(2,:), 'k-.'), grid
title('Controlled Output and Setpoint Profiles (Deviation)')
ylabel('y_i(k)  (deg C)')
legend('y_1(k)', 'r_1(k)', 'y_2(k)', 'r_2(k)', ...
       'Orientation','horizontal','Location','Best')

subplot(212)
plot(kT, Yk(1,:),    'b',  ...
     kT, Rk_f(1,:), 'b-.', ...
     kT, Yk(2,:),    'k',  ...
     kT, Rk_f(2,:), 'k-.'), grid
xlabel('Sampling Instant (k)')
ylabel('Y_i(k)  (deg C)')
legend('T_1(k)', 'r_1(k)', 'T_2(k)', 'r_2(k)', ...
       'Orientation','horizontal','Location','Best')
title('Controlled Output and Setpoint Profiles')

%% ── 15.  Figure 2 – Manipulated input profiles ───────────────────────────
%  Structure mirrors the LQOC Figure 2.
%  In the LQOC script the "target input" us_k is the feed-forward steady-
%  state target.  For PI there is no explicit target input computation, so
%  we display the PI output uk alongside the hard limits (same role).
figure
subplot(211)
stairs(kT, uk(1,:)', 'b',  'LineWidth',1), hold on, grid
stairs(kT, uk(2,:)', 'k',  'LineWidth',1)
stairs(kT, u_limit(:,1),   '-.r','LineWidth',1)   % lower limit
stairs(kT, u_limit(:,2),   '-.r','LineWidth',1)   % upper limit
hold off
ylabel('u_i(k)  (perturbation, %)')
title('Manipulated Input Profiles – Deviation Form')
legend('u_1(k)', 'u_2(k)', 'Limits', ...
       'Orientation','horizontal','Location','Best')

subplot(212)
stairs(kT, Uk(1,:)', 'b', 'LineWidth',1), hold on, grid
stairs(kT, Uk(2,:)', 'k', 'LineWidth',1)
stairs(kT, (u_limit(:,1) + Us(1)*ones(N,1)), '-.r','LineWidth',1)
stairs(kT, (u_limit(:,2) + Us(1)*ones(N,1)), '-.r','LineWidth',1)
hold off
ylabel('U_i(k)  (%  heater power)')
xlabel('Sampling Instant (k)')
title('Manipulated Input Profiles')
legend('U_1(k)', 'U_2(k)', 'Limits', ...
       'Orientation','horizontal','Location','Best')

%% ── 16.  Figure 3 – Model-plant mismatch signal ─────────────────────────
%  In the LQOC script this is LQR.ek_f – the filtered innovation signal.
%  For PI we use ek_f: the low-pass filtered output error, which plays the
%  same diagnostic role (reveals steady-state offset and noise level).
figure
stairs(kT, ek_f', 'LineWidth',1), grid
xlabel('Sampling Instant (k)')
ylabel('e_{f,i}(k)  (deg C)')
title('Filtered Output Error / Model-Plant Mismatch Signal')
legend('e_{f,1}(k)', 'e_{f,2}(k)', ...
       'Orientation','horizontal','Location','Best')

%% ── 17.  Figure 4 – System state profiles ───────────────────────────────
%  Mirrors LQOC Figure 4: absolute state trajectories.
figure
plot(kT, Xk(1,:), 'b', kT, Xk(2,:), 'k'), grid
xlabel('Sampling Instant (k)')
ylabel('x_i(k)  (deg C)')
legend('T_1(k)', 'T_2(k)', 'Orientation','horizontal','Location','Best')
title('System State Profiles')

%% ── 18.  Figure 5 – PI internal signals ─────────────────────────────────
%  Replaces the LQOC "Target states vs Model estimates" plot with PI-
%  specific diagnostics: filtered setpoint and raw vs filtered error.
figure
subplot(211)
plot(kT, rk_f(1,:), 'b-.', kT, yk(1,:), 'b', ...
     kT, rk_f(2,:), 'k-.', kT, yk(2,:), 'k'), grid
ylabel('Signal (deg C)')
title('Filtered Setpoint (dashed) vs Measured Output (solid) – Deviation')
legend('r_{f,1}(k)','y_1(k)','r_{f,2}(k)','y_2(k)', ...
       'Orientation','horizontal','Location','Best')

subplot(212)
plot(kT, ek(1,:),   'b--', ...
     kT, ek_f(1,:), 'b',   ...
     kT, ek(2,:),   'k--', ...
     kT, ek_f(2,:), 'k'), grid
xlabel('Sampling Instant (k)')
ylabel('Error (deg C)')
title('Raw Error (dashed) vs Filtered Error (solid) per PI Loop')
legend('e_1(k)','e_{f,1}(k)','e_2(k)','e_{f,2}(k)', ...
       'Orientation','horizontal','Location','Best')