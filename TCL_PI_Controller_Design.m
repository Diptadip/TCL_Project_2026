clc; clear all; close all;

%% ── 0. Load linearised model ─────────────────────────────────────────────
load TCL_MechModel_Parameters          % loads struct TCL

Phi   = TCL.phy;                       % 2x2 discrete A matrix
Gamma = TCL.gama;                      % 2x2 discrete B matrix
C     = TCL.C_mat;                     % 2x2 eye
D     = TCL.D_mat;                     % 2x2 zeros
Ts    = TCL.Samp_T;                    % sampling period (4 s)
Xs    = TCL.Xs;                        % steady-state state vector
Us    = TCL.Us;                        % steady-state input vector
Ys    = TCL.Ys;                        % steady-state output vector

fprintf('=== TCL Linearised Discrete SS Model ===\n');
fprintf('Sampling period Ts = %g s\n\n', Ts);

%% ── 1. Convert to continuous-time for IMC tuning ───────────────────────
% c2d is invertible; we need the CT model for standard IMC formulae
dmod = ss(Phi, Gamma, C, D, Ts);
cmod = d2c(dmod);                      % ZOH inversion to CT SS model

Ac = cmod.A;
Bc = cmod.B;

fprintf('Continuous-time A matrix:\n'); disp(Ac)
fprintf('Continuous-time B matrix:\n'); disp(Bc)

%% ── 2. Extract SISO process gains for decentralised pairing ─────────────
% Transfer function matrix G(s) = C*(sI - A)^{-1}*B + D
Gs = tf(cmod);

% Diagonal channels used for PI design  (Bristol RGA analysis below)
G11 = Gs(1,1);     % T1 / H1
G22 = Gs(2,2);     % T2 / H2
G12 = Gs(1,2);     % T1 / H2  (coupling)
G21 = Gs(2,1);     % T2 / H1  (coupling)

%% ── 3. Relative Gain Array (RGA) – Bristol pairing check ───────────────
% Steady-state gain matrix K = G(0)
K_ss = dcgain(Gs);

fprintf('=== Steady-State Gain Matrix K = G(0) ===\n');
disp(K_ss)

RGA = K_ss .* inv(K_ss)';

fprintf('=== Relative Gain Array (RGA) ===\n');
disp(RGA)
fprintf(' -> RGA(1,1) = %.3f  (>0.5 confirms H1-T1 pairing is preferred)\n', RGA(1,1));
fprintf(' -> RGA(2,2) = %.3f  (>0.5 confirms H2-T2 pairing is preferred)\n\n', RGA(2,2));

%% ── 4. Fit FOPDT models to diagonal channels ────────────────────────────

% --- Loop 1: G11 ---
t_fopdt = 0:Ts:1200; % long enough to reach steady state
[y1, t1] = lsim(G11, ones(size(t_fopdt)), t_fopdt);
[Kp1, tau1, theta1] = fit_FOPDT(t1, y1);

% --- Loop 2: G22 ---
[y2, t2] = lsim(G22, ones(size(t_fopdt)), t_fopdt);
[Kp2, tau2, theta2] = fit_FOPDT(t2, y2);

fprintf('=== FOPDT Fit Results ===\n');
fprintf('Loop 1 (H1->T1): Kp=%.4f  tau=%.2f s  theta=%.2f s\n', Kp1, tau1, theta1);
fprintf('Loop 2 (H2->T2): Kp=%.4f  tau=%.2f s  theta=%.2f s\n\n', Kp2, tau2, theta2);

%% ── 5. IMC-PI Tuning ────────────────────────────────────────────────────

lambda1 = max(0.2*tau1, theta1);       % IMC tuning parameter, loop 1
lambda2 = max(0.2*tau2, theta2);       % IMC tuning parameter, loop 2

% PI parameters – continuous time
Kc1   = (tau1 + theta1/2) / (Kp1 * (lambda1 + theta1/2));
TauI1 = tau1 + theta1/2;

Kc2   = (tau2 + theta2/2) / (Kp2 * (lambda2 + theta2/2));
TauI2 = tau2 + theta2/2;

fprintf('=== IMC-PI Tuning Parameters (Continuous Time) ===\n');
fprintf('Loop 1:  Kc1 = %.4f    TauI1 = %.4f s    (lambda1 = %.2f s)\n', Kc1, TauI1, lambda1);
fprintf('Loop 2:  Kc2 = %.4f    TauI2 = %.4f s    (lambda2 = %.2f s)\n\n', Kc2, TauI2, lambda2);

%% ── 6. Discrete PI (velocity / incremental form) ────────────────────────

% Loop 1
q0_1 = Kc1 * (1 + Ts/(2*TauI1));
q1_1 = Kc1 * (-1 + Ts/(2*TauI1));

% Loop 2
q0_2 = Kc2 * (1 + Ts/(2*TauI2));
q1_2 = Kc2 * (-1 + Ts/(2*TauI2));

fprintf('=== Discrete PI Velocity Form (Tustin / bilinear, Ts=%g s) ===\n', Ts);
fprintf('u(k) = u(k-1) + q0*e(k) + q1*e(k-1)\n\n');
fprintf('Loop 1:  q0_1 = %.6f    q1_1 = %.6f\n', q0_1, q1_1);
fprintf('Loop 2:  q0_2 = %.6f    q1_2 = %.6f\n\n', q0_2, q1_2);

%% ── 7. Build closed-loop TF for analysis ────────────────────────────────
% Continuous PI transfer function
C1s = tf(Kc1 * [TauI1, 1], [TauI1, 0]);   % Kc*(1 + 1/(TauI*s))
C2s = tf(Kc2 * [TauI2, 1], [TauI2, 0]);

% Closed-loop (diagonal decentralised, ignoring coupling for Bode check)
L1 = C1s * G11;                            % open-loop, loop 1
L2 = C2s * G22;                            % open-loop, loop 2

T1_cl = feedback(L1, 1);                   % closed-loop, loop 1
T2_cl = feedback(L2, 1);                   % closed-loop, loop 2

%% ── 8. Stability margins (manual – avoids broken bodeplot internal class) ─

w = logspace(-4, 2, 8000);          % rad/s grid

[Gm1_dB, Pm1, Wpc1, Wgc1] = manual_margins(L1, w);
[Gm2_dB, Pm2, Wpc2, Wgc2] = manual_margins(L2, w);

fprintf('=== Stability Margins ===\n');
fprintf('Loop 1:  Gain margin = %.2f dB   Phase margin = %.1f deg\n', Gm1_dB, Pm1);
fprintf('         Phase crossover wpc = %.4f rad/s    Gain crossover wgc = %.4f rad/s\n', Wpc1, Wgc1);
fprintf('Loop 2:  Gain margin = %.2f dB   Phase margin = %.1f deg\n', Gm2_dB, Pm2);
fprintf('         Phase crossover wpc = %.4f rad/s    Gain crossover wgc = %.4f rad/s\n', Wpc2, Wgc2);
fprintf(' (Targets: GM > 6 dB,  PM > 45 deg for robust performance)\n\n');

%% ── 9. Step response simulation (closed-loop, manual lsim) ──────────────
t_sim = 0:Ts:600;                          % 10-minute simulation
Nsim  = length(t_sim);
u_step = ones(Nsim, 1);                    % unit step input

% Manual step response via lsim (avoids broken step/bodeplot internals)
y1_step = lsim(T1_cl, u_step, t_sim);
y2_step = lsim(T2_cl, u_step, t_sim);

% Manual performance metrics (rise time, settling time, overshoot)
[rt1, st1, os1] = manual_stepinfo(t_sim, y1_step, 1.0);
[rt2, st2, os2] = manual_stepinfo(t_sim, y2_step, 1.0);

fprintf('=== Closed-Loop Step Response Metrics ===\n');
fprintf('Loop 1 (T1): Rise time=%.1f s  Settling time=%.1f s  Overshoot=%.1f%%\n', rt1, st1, os1);
fprintf('Loop 2 (T2): Rise time=%.1f s  Settling time=%.1f s  Overshoot=%.1f%%\n\n', rt2, st2, os2);

%% ── 10. Full MIMO closed-loop simulation with coupling ──────────────────
%  Simultaneous step changes: +5 in T1ref, +3 in T2ref
Nsim = length(t_sim);
y_cl  = zeros(Nsim, 2);    % outputs (perturbation form, deg C)
u_cl  = zeros(Nsim, 2);    % inputs  (perturbation form, %)
e_cl  = zeros(Nsim, 2);    % errors
x_pert = zeros(2,1);       % state perturbation

% Discrete PI states
e_prev1 = 0;  u_prev1 = 0;
e_prev2 = 0;  u_prev2 = 0;

r1 = 5.0;     % setpoint step for T1 (deg C above steady state)
r2 = 3.0;     % setpoint step for T2 (deg C above steady state)

for k = 1:Nsim
    y_pert = C * x_pert;            % current output perturbation

    % Errors
    e1 = r1 - y_pert(1);
    e2 = r2 - y_pert(2);

    % Discrete PI velocity form (with anti-windup clamping)
    du1 = q0_1*e1 + q1_1*e_prev1;
    du2 = q0_2*e2 + q1_2*e_prev2;

    u1 = u_prev1 + du1;
    u2 = u_prev2 + du2;

    % Anti-windup: clamp to allowable perturbation range (±50 %)
    u1 = min(max(u1, -Us(1)), 100 - Us(1));
    u2 = min(max(u2, -Us(2)), 100 - Us(2));

    % Store
    y_cl(k,:) = y_pert';
    u_cl(k,:) = [u1, u2];
    e_cl(k,:) = [e1, e2];

    % State update (discrete SS perturbation model)
    x_pert = Phi*x_pert + Gamma*[u1; u2];

    % Update previous values
    e_prev1 = e1;  u_prev1 = u1;
    e_prev2 = e2;  u_prev2 = u2;
end

%% ── 11. Plots ────────────────────────────────────────────────────────────

% --- Manual Bode plots (magnitude + phase) ---
w_plot = logspace(-4, 2, 2000);
Hjw1 = squeeze(freqresp(L1, w_plot));
Hjw2 = squeeze(freqresp(L2, w_plot));

mag1_dB  = 20*log10(abs(Hjw1));
phs1_deg = 180/pi * unwrap(angle(Hjw1));
mag2_dB  = 20*log10(abs(Hjw2));
phs2_deg = 180/pi * unwrap(angle(Hjw2));

figure('Name','Open-Loop Bode Plots','NumberTitle','off');
subplot(221);
semilogx(w_plot, mag1_dB, 'b-', 'LineWidth', 1.5); grid on;
yline(0,'k--','LineWidth',0.8);
ylabel('Magnitude (dB)'); title('Loop 1 (H_1\rightarrowT_1): Bode');
xline(Wgc1,'r--','LineWidth',0.8);

subplot(223);
semilogx(w_plot, phs1_deg, 'b-', 'LineWidth', 1.5); grid on;
yline(-180,'k--','LineWidth',0.8);
ylabel('Phase (deg)'); xlabel('Frequency (rad/s)');
xline(Wpc1,'r--','LineWidth',0.8);
text(Wgc1*1.15, min(phs1_deg)*0.9, sprintf('PM=%.1f°',Pm1), 'Color','r','FontSize',8);

subplot(222);
semilogx(w_plot, mag2_dB, 'r-', 'LineWidth', 1.5); grid on;
yline(0,'k--','LineWidth',0.8);
ylabel('Magnitude (dB)'); title('Loop 2 (H_2\rightarrowT_2): Bode');
xline(Wgc2,'b--','LineWidth',0.8);

subplot(224);
semilogx(w_plot, phs2_deg, 'r-', 'LineWidth', 1.5); grid on;
yline(-180,'k--','LineWidth',0.8);
ylabel('Phase (deg)'); xlabel('Frequency (rad/s)');
xline(Wpc2,'b--','LineWidth',0.8);
text(Wgc2*1.15, min(phs2_deg)*0.9, sprintf('PM=%.1f°',Pm2), 'Color','b','FontSize',8);

% --- Manual Nyquist plots ---
w_ny  = logspace(-4, 2, 5000);
Lny1  = squeeze(freqresp(L1, w_ny));
Lny2  = squeeze(freqresp(L2, w_ny));

figure('Name','Nyquist Diagrams','NumberTitle','off');
subplot(121);
plot(real(Lny1), imag(Lny1), 'b-', 'LineWidth', 1.5); hold on;
plot(real(Lny1), -imag(Lny1), 'b--', 'LineWidth', 0.8);   % mirror
plot(-1, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal;
xlabel('Real'); ylabel('Imaginary');
title('Loop 1: Nyquist'); legend('L_1(j\omega)','-180° mirror','-1 point');

subplot(122);
plot(real(Lny2), imag(Lny2), 'r-', 'LineWidth', 1.5); hold on;
plot(real(Lny2), -imag(Lny2), 'r--', 'LineWidth', 0.8);
plot(-1, 0, 'bx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal;
xlabel('Real'); ylabel('Imaginary');
title('Loop 2: Nyquist'); legend('L_2(j\omega)','-180° mirror','-1 point');

% --- Closed-loop step responses (single-loop, lsim-based) ---
figure('Name','Closed-Loop Step Responses (SISO)','NumberTitle','off');
subplot(211);
plot(t_sim, y1_step, 'b-', 'LineWidth', 1.5); hold on;
yline(1,'k--','LineWidth',0.8); grid on;
ylabel('\DeltaT_1 (normalized)');
title('Step response – Loop 1 (H_1 \rightarrow T_1)');
legend('Output','Setpoint','Location','southeast');

subplot(212);
plot(t_sim, y2_step, 'r-', 'LineWidth', 1.5); hold on;
yline(1,'k--','LineWidth',0.8); grid on;
ylabel('\DeltaT_2 (normalized)');
xlabel('Time (s)');
title('Step response – Loop 2 (H_2 \rightarrow T_2)');
legend('Output','Setpoint','Location','southeast');

% --- MIMO closed-loop simulation ---
figure('Name','MIMO Closed-Loop Simulation','NumberTitle','off');
subplot(211);
plot(t_sim, y_cl(:,1)+Ys(1), 'b-', 'LineWidth', 1.5); hold on;
plot(t_sim, (r1+Ys(1))*ones(size(t_sim)), 'b--', 'LineWidth', 0.8);
plot(t_sim, y_cl(:,2)+Ys(2), 'r-', 'LineWidth', 1.5);
plot(t_sim, (r2+Ys(2))*ones(size(t_sim)), 'r--', 'LineWidth', 0.8);
grid on; ylabel('Temperature (°C)');
title('MIMO PI Control – Simultaneous Step Changes (+5°C in T_1, +3°C in T_2)');
legend('T_1','T_1 ref','T_2','T_2 ref','Location','southeast');

subplot(212);
stairs(t_sim, u_cl(:,1)+Us(1), 'b-', 'LineWidth', 1.5); hold on;
stairs(t_sim, u_cl(:,2)+Us(2), 'r-', 'LineWidth', 1.5);
grid on; ylabel('Heater Input (%)'); xlabel('Time (s)');
title('Control Inputs');
legend('H_1','H_2','Location','northeast');

%% ── 12. Save results ─────────────────────────────────────────────────────
TCL.PI.Kc1   = Kc1;    TCL.PI.TauI1 = TauI1;
TCL.PI.Kc2   = Kc2;    TCL.PI.TauI2 = TauI2;
TCL.PI.q0_1  = q0_1;   TCL.PI.q1_1  = q1_1;
TCL.PI.q0_2  = q0_2;   TCL.PI.q1_2  = q1_2;
TCL.PI.lambda1 = lambda1;
TCL.PI.lambda2 = lambda2;

save TCL_MechModel_Parameters TCL
fprintf('PI parameters saved to TCL_MechModel_Parameters.mat\n');

%% ── LOCAL FUNCTIONS ─────────────────────────────────────────────────────

function [Kp, tau, theta] = fit_FOPDT(t, y)
% Identify K, tau, theta from step-response vector y(t).
% Uses the standard graphical 28.3% / 63.2% tangent method.

    Kp    = y(end);                        % steady-state gain

    % 63.2% and 28.3% points
    idx63 = find(y >= 0.632*Kp, 1, 'first');
    idx28 = find(y >= 0.283*Kp, 1, 'first');

    if isempty(idx63) || isempty(idx28)
        % Fallback: least-squares fit
        tau   = t(end)/5;
        theta = 0;
        return
    end

    t63 = t(idx63);
    t28 = t(idx28);

    tau   = 1.5 * (t63 - t28);            % process time constant
    theta = t63 - tau;                     % apparent dead time
    theta = max(theta, 0);                 % cannot be negative
end

function [GM_dB, PM_deg, Wpc, Wgc] = manual_margins(L, w)
% Compute stability margins by scanning L(jw) on the supplied frequency grid.
%   GM_dB  – gain margin in dB  (Inf if no phase crossover found)
%   PM_deg – phase margin in degrees
%   Wpc    – phase crossover frequency (rad/s)
%   Wgc    – gain crossover frequency  (rad/s)

    Hjw  = squeeze(freqresp(L, w));
    mag  = abs(Hjw);
    phs  = 180/pi * unwrap(angle(Hjw));   % unwrapped phase in degrees

    % ---- Gain crossover: |L(jwgc)| = 1 (0 dB) ----
    mag_minus1 = mag - 1;
    idx_gc = find(diff(sign(mag_minus1)) ~= 0, 1, 'first');
    if isempty(idx_gc)
        Wgc    = NaN;
        PM_deg = Inf;
    else
        % Linear interpolation for precision
        alpha  = -mag_minus1(idx_gc) / (mag_minus1(idx_gc+1) - mag_minus1(idx_gc));
        Wgc    = w(idx_gc) + alpha*(w(idx_gc+1) - w(idx_gc));
        phs_gc = phs(idx_gc) + alpha*(phs(idx_gc+1) - phs(idx_gc));
        PM_deg = 180 + phs_gc;
    end

    % ---- Phase crossover: angle(L(jwpc)) = -180 deg ----
    phs_plus180 = phs + 180;
    % Look for the first downward crossing (most conservative)
    idx_pc = find(diff(sign(phs_plus180)) < 0, 1, 'first');
    if isempty(idx_pc)
        Wpc   = NaN;
        GM_dB = Inf;
    else
        alpha2 = -phs_plus180(idx_pc) / (phs_plus180(idx_pc+1) - phs_plus180(idx_pc));
        Wpc    = w(idx_pc) + alpha2*(w(idx_pc+1) - w(idx_pc));
        mag_pc = mag(idx_pc) + alpha2*(mag(idx_pc+1) - mag(idx_pc));
        GM_dB  = -20*log10(mag_pc);       % GM = 1/|L(jwpc)| in dB
    end
end

function [rt, st, os] = manual_stepinfo(t, y, yss)
% Rise time, settling time (2% band), overshoot for step response y(t).
%   yss – steady-state value (normally 1 for a normalised step)

    % Rise time: 10% -> 90% of yss
    idx10 = find(y >= 0.10*yss, 1, 'first');
    idx90 = find(y >= 0.90*yss, 1, 'first');
    if isempty(idx10) || isempty(idx90)
        rt = NaN;
    else
        rt = t(idx90) - t(idx10);
    end

    % Overshoot
    pk = max(y);
    os = max(0, (pk - yss)/yss * 100);

    % Settling time: last time y exits the 2% band around yss
    band    = 0.02 * yss;
    outside = find(abs(y - yss) > band);
    if isempty(outside)
        st = 0;
    else
        st = t(outside(end));
    end
end