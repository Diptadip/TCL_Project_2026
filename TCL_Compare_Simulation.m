% =========================================================================
%  TCL_Compare_Simulation.m
%
%  PURPOSE
%  -------
%  Runs BOTH the PI servo simulation AND the MPC black-box simulation
%  (noise-free only) on the NONLINEAR / linear plant models respectively,
%  then produces a graphical comparison and a performance-index table.
%
%  REQUIREMENTS BEFORE RUNNING
%  ---------------------------
%    1. TCL_MechModel_Parameters.mat   (TCL struct with TCL.PI fields)
%    2. TCL_BlackBox_OE_SS.mat         (idmod struct)
%    Both files must be on the MATLAB path.
%
%  COMPARISON PLOTS
%  ----------------
%    Fig 1  –  Yi(k) vs k and Ri(k) vs k  for PI and MPC  (i=1,2)
%    Fig 2  –  Ui(k) vs k  (stairs)  for PI and MPC  (i=1,2)
%    Fig 3  –  ef,i(k) vs k  (stairs)  for PI and MPC  (i=1,2)
%    Fig 4  –  xs,i(k) and us,i(k) vs k  for MPC only  (i=1,2)
%
%  PERFORMANCE INDICES  (printed as a formatted table)
%  ---------------------------------------------------
%    SSE_i   = sum_{k=1}^{Ns} [Yi(k) - Ri(k)]^2          i=1,2
%    SSMV_i  = sum_{k=1}^{Ns} [Ui(k) - Us_i]^2            i=1,2
%    SSdelMV_i= sum_{k=1}^{Ns} [Ui(k) - Ui(k-1)]^2        i=1,2
%
%  OUTPUT
%    TCL_Compare_SimResults.mat
% =========================================================================

clear all
close all
clc

global TCL

%% ── 0. Load models ───────────────────────────────────────────────────────
load TCL_MechModel_Parameters.mat    % TCL struct
load TCL_BlackBox_OE_SS.mat          % idmod struct

if ~isfield(TCL,'PI')
    error('TCL.PI not found. Run TCL_PI_Controller_Design.m first.');
end

% Operating point
Ys = TCL.Ys(:);   % 2x1  deg C
Us = TCL.Us(:);   % 2x1  %
Ts = TCL.Samp_T;  % 4 s

% Black-box model
A_bb = idmod.phy;
B_bb = idmod.gama;
C_bb = idmod.C_mat;
D_bb = idmod.D_mat;

[n_op, n_st_bb] = size(C_bb);
[~,    n_ip    ] = size(B_bb);
[~,    n_st_pi ] = size(TCL.gama);

%% ── 1. Common experiment timing ──────────────────────────────────────────
N        = 550;
k_warmup = 30;
k_settle = 100;
k_step1  = 101;
k_step2  = 301;

delta_r_step1 = [ +6; -4 ];
delta_r_step2 = [ -5; +7 ];

kT    = (0:N-1)';
t_sec = kT * Ts;

u_H = [100; 100] - Us;
u_L = [  0;   0] - Us;

%% ── 2. Filter parameters (same for both controllers) ────────────────────
beta_r = 0.95;
alfa_e = 0.95;
phy_r  = beta_r * eye(n_op);
phy_e  = alfa_e * eye(n_op);

%% ════════════════════════════════════════════════════════════════════════
%                        PI  SIMULATION
%  Uses NONLINEAR plant (ode45 + TCL_Dynamics) exactly as in
%  TCL_PI_ServoDesign_Simulation.m  (Plant_Sim = 0)
% ════════════════════════════════════════════════════════════════════════

fprintf('Running PI closed-loop simulation (nonlinear plant) ...\n');

q0_1 = TCL.PI.q0_1;  q1_1 = TCL.PI.q1_1;
q0_2 = TCL.PI.q0_2;  q1_2 = TCL.PI.q1_2;

% Allocate
xk_pi   = zeros(n_st_pi, N);
uk_pi   = zeros(n_ip,    N);
yk_pi   = zeros(n_op,    N);
rk_f_pi = zeros(n_op,    N);
ek_f_pi = zeros(n_op,    N);

% Noise (off for fair comparison)
vk_pi = zeros(n_op, N);

e_prev_pi = zeros(n_ip, 1);
u_prev_pi = zeros(n_ip, 1);
rk_f_now_pi = zeros(n_op, 1);
ek_f_now_pi = zeros(n_op, 1);

yk_pi(:,1) = TCL.C_mat * xk_pi(:,1);

for k = 1:N-1

    % Setpoint schedule
    if k < k_step1
        rk = zeros(n_op,1);
    elseif k < k_step2
        rk = delta_r_step1;
    else
        rk = delta_r_step2;
    end

    % Setpoint filter
    if k == 1
        rk_f_now_pi = (eye(n_op) - phy_r) * rk;
    else
        rk_f_now_pi = phy_r * rk_f_pi(:,k-1) + (eye(n_op) - phy_r) * rk;
    end
    rk_f_pi(:,k) = rk_f_now_pi;

    % Error and innovation filter
    raw_ek = rk_f_now_pi - yk_pi(:,k);
    if k == 1
        ek_f_now_pi = (eye(n_op) - phy_e) * raw_ek;
    else
        ek_f_now_pi = phy_e * ek_f_pi(:,k-1) + (eye(n_op) - phy_e) * raw_ek;
    end
    ek_f_pi(:,k) = ek_f_now_pi;

    % Warm-up
    if k <= k_warmup
        u1 = 0;  u2 = 0;
        e_prev_pi = zeros(n_ip,1);
        u_prev_pi = zeros(n_ip,1);
    else
        % PI velocity law (Tustin)
        du1 = q0_1 * ek_f_now_pi(1) + q1_1 * e_prev_pi(1);
        du2 = q0_2 * ek_f_now_pi(2) + q1_2 * e_prev_pi(2);
        u1  = u_prev_pi(1) + du1;
        u2  = u_prev_pi(2) + du2;

        % Anti-windup
        u1 = min(max(u1, u_L(1)), u_H(1));
        u2 = min(max(u2, u_L(2)), u_H(2));

        e_prev_pi = ek_f_now_pi;
        u_prev_pi = [u1; u2];
    end
    uk_pi(:,k) = [u1; u2];

    % Plant: nonlinear ODE
    TCL.Uk = TCL.Us + uk_pi(:,k);
    [~, Xt] = ode45('TCL_Dynamics', [0 Ts], TCL.Xs + xk_pi(:,k));
    xk_pi(:,k+1) = Xt(end,:)' - TCL.Xs;

    % Output
    yk_pi(:,k+1) = TCL.C_mat * xk_pi(:,k+1) + vk_pi(:,k+1);
end

% Pad last sample
uk_pi(:,N)   = uk_pi(:,N-1);
ek_f_pi(:,N) = ek_f_pi(:,N-1);
rk_f_pi(:,N) = rk_f_pi(:,N-1);

% Absolute arrays
Yk_pi   = yk_pi   + repmat(Ys, 1, N);
Uk_pi   = uk_pi   + repmat(Us, 1, N);
Rk_f_pi = rk_f_pi + repmat(Ys, 1, N);

fprintf('PI simulation done.\n\n');

%% ════════════════════════════════════════════════════════════════════════
%                        MPC  SIMULATION  (noise-free)
%  Uses black-box linear model as plant (same as TCL_MPC_BB_Simulation.m,
%  Noise_ON = 0)
% ════════════════════════════════════════════════════════════════════════

fprintf('Running MPC closed-loop simulation (black-box linear plant, noise-free) ...\n');

%% ── Kalman filter design ─────────────────────────────────────────────────
Q_kf = 1e-3 * eye(n_st_bb);
R_kf = TCL.R;
[~, ~, L_kf] = dare(A_bb', C_bb', Q_kf, R_kf);
L_kf = L_kf';

%% ── Steady-state target gain ─────────────────────────────────────────────
M_ss     = [(eye(n_st_bb)-A_bb), -B_bb; C_bb, zeros(n_op, n_ip)];
M_ss_inv = pinv(M_ss);

%% ── MPC tuning ───────────────────────────────────────────────────────────
Np    = 15;
Nc    = 5;
Wx    = diag([15,  15 ]);
Wdelu = diag([1.0, 1.0]);

n_aug = n_st_bb + n_ip;
A_aug = [A_bb, B_bb; zeros(n_ip, n_st_bb), eye(n_ip)];
B_aug = [B_bb; eye(n_ip)];
C_aug = [C_bb, zeros(n_op, n_ip)];

Psi   = zeros(Np*n_op, n_aug);
Theta = zeros(Np*n_op, Nc*n_ip);

A_aug_pow = eye(n_aug);
for i = 1:Np
    A_aug_pow      = A_aug_pow * A_aug;
    rows_i         = (i-1)*n_op + 1 : i*n_op;
    Psi(rows_i,:)  = C_aug * A_aug_pow;
    for j = 1:min(i,Nc)
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
K_mpc     = K_mpc(1:n_ip, :);

%% ── MPC simulation loop ──────────────────────────────────────────────────
xk_hat_mpc  = zeros(n_st_bb, N);
uk_mpc      = zeros(n_ip,    N);
yk_mpc      = zeros(n_op,    N);
rk_f_mpc    = zeros(n_op,    N);
ek_f_mpc    = zeros(n_op,    N);
xs_mpc      = zeros(n_st_bb, N);
us_mpc      = zeros(n_ip,    N);
xk_true_mpc = zeros(n_st_bb, N);

u_prev_mpc = zeros(n_ip, 1);

for k = 1:N

    % Setpoint schedule
    if k < k_step1
        rk = zeros(n_op,1);
    elseif k < k_step2
        rk = delta_r_step1;
    else
        rk = delta_r_step2;
    end

    % Setpoint filter
    if k == 1
        rk_f_mpc(:,k) = (eye(n_op) - phy_r) * rk;
    else
        rk_f_mpc(:,k) = phy_r * rk_f_mpc(:,k-1) + (eye(n_op) - phy_r) * rk;
    end

    % Error and innovation filter
    raw_ek = rk_f_mpc(:,k) - yk_mpc(:,k);
    if k == 1
        ek_f_mpc(:,k) = (eye(n_op) - phy_e) * raw_ek;
    else
        ek_f_mpc(:,k) = phy_e * ek_f_mpc(:,k-1) + (eye(n_op) - phy_e) * raw_ek;
    end

    % Steady-state targets
    rhs           = [zeros(n_st_bb,1); rk_f_mpc(:,k)];
    xus           = M_ss_inv * rhs;
    xs_mpc(:,k)   = xus(1:n_st_bb);
    us_mpc(:,k)   = xus(n_st_bb+1:end);

    % Warm-up
    if k <= k_warmup
        uk_mpc(:,k)        = zeros(n_ip,1);
        u_prev_mpc         = zeros(n_ip,1);
        xk_hat_mpc(:,k)    = zeros(n_st_bb,1);
    else
        x_aug  = [xk_hat_mpc(:,k); u_prev_mpc];
        R_vec  = repmat(rk_f_mpc(:,k), Np, 1);
        du_opt = K_mpc * (R_vec - Psi * x_aug);
        u_new  = u_prev_mpc + du_opt;
        u_new  = min(max(u_new, u_L), u_H);
        uk_mpc(:,k) = u_new;
        u_prev_mpc  = u_new;
    end

    % Plant update (linear black-box, noise-free)
    if k < N
        xk_true_mpc(:,k+1) = A_bb * xk_true_mpc(:,k) + B_bb * uk_mpc(:,k);
        yk_mpc(:,k+1)      = C_bb * xk_true_mpc(:,k+1);

        % Kalman update
        y_pred            = C_bb * xk_hat_mpc(:,k);
        inno              = yk_mpc(:,k) - y_pred;
        xk_hat_mpc(:,k+1) = A_bb * xk_hat_mpc(:,k) + B_bb * uk_mpc(:,k) + L_kf * inno;
    end
end

% Absolute arrays
Yk_mpc   = yk_mpc   + repmat(Ys, 1, N);
Uk_mpc   = uk_mpc   + repmat(Us, 1, N);
Rk_f_mpc = rk_f_mpc + repmat(Ys, 1, N);
Xs_mpc   = xs_mpc   + repmat(TCL.Xs, 1, N);
Us_mpc   = us_mpc   + repmat(Us,     1, N);

fprintf('MPC simulation done.\n\n');

%% ════════════════════════════════════════════════════════════════════════
%                     PERFORMANCE INDICES
% ════════════════════════════════════════════════════════════════════════
% Use full simulation window Ns = N

Ns = N;

% --- SSE ---
SSE_pi_1  = sum((Yk_pi(1,1:Ns)   - Rk_f_pi(1,1:Ns)).^2);
SSE_pi_2  = sum((Yk_pi(2,1:Ns)   - Rk_f_pi(2,1:Ns)).^2);
SSE_mpc_1 = sum((Yk_mpc(1,1:Ns)  - Rk_f_mpc(1,1:Ns)).^2);
SSE_mpc_2 = sum((Yk_mpc(2,1:Ns)  - Rk_f_mpc(2,1:Ns)).^2);

% --- SSMV ---
SSMV_pi_1  = sum((Uk_pi(1,1:Ns)  - Us(1)).^2);
SSMV_pi_2  = sum((Uk_pi(2,1:Ns)  - Us(2)).^2);
SSMV_mpc_1 = sum((Uk_mpc(1,1:Ns) - Us(1)).^2);
SSMV_mpc_2 = sum((Uk_mpc(2,1:Ns) - Us(2)).^2);

% --- SSdelMV ---
dUk_pi_1  = diff(Uk_pi(1,1:Ns));
dUk_pi_2  = diff(Uk_pi(2,1:Ns));
dUk_mpc_1 = diff(Uk_mpc(1,1:Ns));
dUk_mpc_2 = diff(Uk_mpc(2,1:Ns));

SSdelMV_pi_1  = sum(dUk_pi_1.^2);
SSdelMV_pi_2  = sum(dUk_pi_2.^2);
SSdelMV_mpc_1 = sum(dUk_mpc_1.^2);
SSdelMV_mpc_2 = sum(dUk_mpc_2.^2);

%% ── Print performance table ───────────────────────────────────────────────
fprintf('=================================================================\n');
fprintf('           SIMULATION PERFORMANCE INDEX COMPARISON\n');
fprintf('           (Ns = %d samples,  %.1f min)\n', Ns, Ns*Ts/60);
fprintf('=================================================================\n');
fprintf('  Index          |   PI (Loop 1) |  PI (Loop 2) | MPC (Loop 1) | MPC (Loop 2)\n');
fprintf('-----------------+---------------+--------------+--------------+-------------\n');
fprintf('  SSE            |  %11.2f  |  %10.2f  |  %10.2f  | %10.2f\n', ...
    SSE_pi_1,  SSE_pi_2,  SSE_mpc_1,  SSE_mpc_2);
fprintf('  SSMV           |  %11.2f  |  %10.2f  |  %10.2f  | %10.2f\n', ...
    SSMV_pi_1, SSMV_pi_2, SSMV_mpc_1, SSMV_mpc_2);
fprintf('  SSdeltaMV      |  %11.2f  |  %10.2f  |  %10.2f  | %10.2f\n', ...
    SSdelMV_pi_1, SSdelMV_pi_2, SSdelMV_mpc_1, SSdelMV_mpc_2);
fprintf('=================================================================\n\n');

%% ════════════════════════════════════════════════════════════════════════
%                          COMPARISON PLOTS
% ════════════════════════════════════════════════════════════════════════

step_locs = [k_step1-1, k_step2-1];

%% ── Figure 1: Yi(k) vs k and Ri(k) vs k ─────────────────────────────────
figure('Name','Comparison – Output Profiles','NumberTitle','off','Position',[50 50 1000 700])

subplot(2,1,1)
plot(kT, Yk_pi(1,:),   'b-',  'LineWidth',1.8), hold on
plot(kT, Yk_mpc(1,:),  'r-',  'LineWidth',1.8)
plot(kT, Rk_f_pi(1,:), 'k--', 'LineWidth',1.2)
xline(step_locs(1),'k:','Step 1','LabelVerticalAlignment','bottom','FontSize',8)
xline(step_locs(2),'k:','Step 2','LabelVerticalAlignment','bottom','FontSize',8)
grid on
ylabel('T_1  (deg C)','FontSize',11)
title('Output Y_1(k) and Reference R_1(k) – Simulation (Noise-Free)','FontSize',12)
legend('Y_1^{PI}(k)','Y_1^{MPC}(k)','R_1(k)','Location','Best','FontSize',9)

subplot(2,1,2)
plot(kT, Yk_pi(2,:),   'b-',  'LineWidth',1.8), hold on
plot(kT, Yk_mpc(2,:),  'r-',  'LineWidth',1.8)
plot(kT, Rk_f_pi(2,:), 'k--', 'LineWidth',1.2)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('T_2  (deg C)','FontSize',11)
xlabel('Sample  k','FontSize',11)
title('Output Y_2(k) and Reference R_2(k)','FontSize',12)
legend('Y_2^{PI}(k)','Y_2^{MPC}(k)','R_2(k)','Location','Best','FontSize',9)

%% ── Figure 2: Ui(k) vs k ─────────────────────────────────────────────────
figure('Name','Comparison – Heater Inputs','NumberTitle','off','Position',[100 50 1000 700])

subplot(2,1,1)
stairs(kT, Uk_pi(1,:)',  'b-', 'LineWidth',1.8), hold on
stairs(kT, Uk_mpc(1,:)', 'r-', 'LineWidth',1.8)
yline(5,  'k--', 'LineWidth',0.9)
yline(80, 'k--', 'LineWidth',0.9)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylim([0 105])
ylabel('U_1(k)  (% heater)','FontSize',11)
title('Heater Input U_1(k) – Simulation','FontSize',12)
legend('U_1^{PI}(k)','U_1^{MPC}(k)','Safety bounds','Location','Best','FontSize',9)

subplot(2,1,2)
stairs(kT, Uk_pi(2,:)',  'b-', 'LineWidth',1.8), hold on
stairs(kT, Uk_mpc(2,:)', 'r-', 'LineWidth',1.8)
yline(5,  'k--', 'LineWidth',0.9)
yline(80, 'k--', 'LineWidth',0.9)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylim([0 105])
ylabel('U_2(k)  (% heater)','FontSize',11)
xlabel('Sample  k','FontSize',11)
title('Heater Input U_2(k)','FontSize',12)
legend('U_2^{PI}(k)','U_2^{MPC}(k)','Safety bounds','Location','Best','FontSize',9)

%% ── Figure 3: ef(k) vs k ─────────────────────────────────────────────────
figure('Name','Comparison – Filtered Error','NumberTitle','off','Position',[150 50 1000 700])

subplot(2,1,1)
stairs(kT, ek_f_pi(1,:)',  'b-', 'LineWidth',1.8), hold on
stairs(kT, ek_f_mpc(1,:)', 'r-', 'LineWidth',1.8)
yline(0,'k--','LineWidth',0.8)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('e_{f,1}(k)  (deg C)','FontSize',11)
title('Filtered Error e_{f,1}(k) – Simulation','FontSize',12)
legend('PI','MPC','Location','Best','FontSize',9)

subplot(2,1,2)
stairs(kT, ek_f_pi(2,:)',  'b-', 'LineWidth',1.8), hold on
stairs(kT, ek_f_mpc(2,:)', 'r-', 'LineWidth',1.8)
yline(0,'k--','LineWidth',0.8)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('e_{f,2}(k)  (deg C)','FontSize',11)
xlabel('Sample  k','FontSize',11)
title('Filtered Error e_{f,2}(k)','FontSize',12)
legend('PI','MPC','Location','Best','FontSize',9)

%% ── Figure 4: xs(k) and us(k) vs k  (MPC only) ──────────────────────────
figure('Name','MPC – SS Targets xs(k) and us(k)','NumberTitle','off','Position',[200 50 1000 700])

subplot(2,1,1)
plot(kT, Xs_mpc(1,:), 'b-', 'LineWidth',1.8), hold on
if n_st_bb > 1
    plot(kT, Xs_mpc(2,:), 'r-', 'LineWidth',1.8)
end
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('x_{s,i}(k)  (deg C)','FontSize',11)
title('MPC SS Target States x_s(k)  [Simulation]','FontSize',12)
leg_xs = arrayfun(@(i) sprintf('x_{s,%d}^{MPC}(k)',i), 1:n_st_bb, 'UniformOutput',false);
legend(leg_xs{:},'Location','Best','FontSize',9)

subplot(2,1,2)
stairs(kT, Us_mpc(1,:)', 'b-', 'LineWidth',1.8), hold on
stairs(kT, Us_mpc(2,:)', 'r-', 'LineWidth',1.8)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('u_{s,i}(k)  (% heater)','FontSize',11)
xlabel('Sample  k','FontSize',11)
title('MPC SS Target Inputs u_s(k)','FontSize',12)
legend('u_{s,1}^{MPC}(k)','u_{s,2}^{MPC}(k)','Location','Best','FontSize',9)

%% ── Figure 5: Performance index bar chart ───────────────────────────────
figure('Name','Performance Index Comparison','NumberTitle','off','Position',[250 50 1000 500])

metrics_loop1 = [SSE_pi_1,     SSE_mpc_1;
                 SSMV_pi_1,    SSMV_mpc_1;
                 SSdelMV_pi_1, SSdelMV_mpc_1];

metrics_loop2 = [SSE_pi_2,     SSE_mpc_2;
                 SSMV_pi_2,    SSMV_mpc_2;
                 SSdelMV_pi_2, SSdelMV_mpc_2];

labels = {'SSE','SSMV','SS\DeltaMV'};
group  = categorical(labels);
group  = reordercats(group, labels);

subplot(1,2,1)
b1 = bar(group, metrics_loop1, 'grouped');
b1(1).FaceColor = [0.2 0.4 0.8];
b1(2).FaceColor = [0.8 0.2 0.2];
grid on
ylabel('Index Value','FontSize',10)
title('Loop 1  (T_1 / H_1)','FontSize',11)
legend('PI','MPC','Location','Best','FontSize',9)

subplot(1,2,2)
b2 = bar(group, metrics_loop2, 'grouped');
b2(1).FaceColor = [0.2 0.4 0.8];
b2(2).FaceColor = [0.8 0.2 0.2];
grid on
ylabel('Index Value','FontSize',10)
title('Loop 2  (T_2 / H_2)','FontSize',11)
legend('PI','MPC','Location','Best','FontSize',9)
sgtitle('Performance Indices – Simulation Comparison','FontSize',13)

%% ── Save ─────────────────────────────────────────────────────────────────
% Pack results
sim_compare.kT          = kT;
sim_compare.t_sec       = t_sec;
sim_compare.Yk_pi       = Yk_pi;
sim_compare.Uk_pi       = Uk_pi;
sim_compare.Rk_f_pi     = Rk_f_pi;
sim_compare.ek_f_pi     = ek_f_pi;
sim_compare.Yk_mpc      = Yk_mpc;
sim_compare.Uk_mpc      = Uk_mpc;
sim_compare.Rk_f_mpc    = Rk_f_mpc;
sim_compare.ek_f_mpc    = ek_f_mpc;
sim_compare.Xs_mpc      = Xs_mpc;
sim_compare.Us_mpc      = Us_mpc;
sim_compare.SSE_pi      = [SSE_pi_1;  SSE_pi_2 ];
sim_compare.SSE_mpc     = [SSE_mpc_1; SSE_mpc_2];
sim_compare.SSMV_pi     = [SSMV_pi_1;  SSMV_pi_2 ];
sim_compare.SSMV_mpc    = [SSMV_mpc_1; SSMV_mpc_2];
sim_compare.SSdelMV_pi  = [SSdelMV_pi_1;  SSdelMV_pi_2 ];
sim_compare.SSdelMV_mpc = [SSdelMV_mpc_1; SSdelMV_mpc_2];

save TCL_Compare_SimResults  sim_compare
fprintf('Results saved to  TCL_Compare_SimResults.mat\n');
fprintf('All done.\n');
