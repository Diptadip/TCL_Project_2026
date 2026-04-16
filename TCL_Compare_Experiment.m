clear all
close all
clc

%% ── 0. Load operating-point info ─────────────────────────────────────────
load TCL_MechModel_Parameters.mat   % TCL struct  (Ys, Us, Samp_T)
Ys = TCL.Ys(:);
Us = TCL.Us(:);
Ts = TCL.Samp_T;   % 4 s

%% ── 1. Load PI experiment results ────────────────────────────────────────
%  Variables expected from TCL_PI_ServoExperiment.m:
%    time, t1s, t2s, h1s, h2s, R1s, R2s, e1s, e2s
load TCL_PI_Servo_ExpResults.mat

% Rename to common convention  (column vectors, absolute values)
kT_pi    = (0:length(t1s)-1)';
Yk_pi    = [t1s(:)'; t2s(:)'];    % 2 × Npi
Uk_pi    = [h1s(:)'; h2s(:)'];
Rk_f_pi  = [R1s(:)'; R2s(:)'];
ek_f_pi  = [e1s(:)'; e2s(:)'];

Npi = length(kT_pi);

fprintf('PI experiment:  %d samples  (%.1f min)\n', Npi, Npi*Ts/60);

%% ── 2. Load MPC experiment results ───────────────────────────────────────
%  Variables expected from TCL_MPC_BB_Experiment.m:
%    kT, t1s, t2s, h1s, h2s, R1s, R2s, e1s, e2s, xs1s, xs2s, us1s, us2s
load TCL_MPC_BB_ExpResults.mat

% Resolve naming conflict: re-assign after PI load
kT_mpc   = kT(:);
Yk_mpc   = [t1s(:)'; t2s(:)'];
Uk_mpc   = [h1s(:)'; h2s(:)'];
Rk_f_mpc = [R1s(:)'; R2s(:)'];
ek_f_mpc = [e1s(:)'; e2s(:)'];

% SS targets (MPC only)
xs_mpc   = [xs1s(:)'; xs2s(:)'];          % perturbation form
us_mpc   = [us1s(:)'; us2s(:)'];          % perturbation form
Xs_mpc   = xs_mpc + repmat(TCL.Xs, 1, length(kT_mpc));   % absolute
Us_mpc   = us_mpc + repmat(Us,     1, length(kT_mpc));    % absolute

Nmpc = length(kT_mpc);

fprintf('MPC experiment: %d samples  (%.1f min)\n', Nmpc, Nmpc*Ts/60);

%% ── 3. Align lengths for fair comparison ────────────────────────────────
%  Use the shorter of the two datasets as the comparison window
Ns = min(Npi, Nmpc);
fprintf('Comparison window  Ns = %d samples  (%.1f min)\n\n', Ns, Ns*Ts/60);

kT_cmp = (0:Ns-1)';

%% ── 4. Step transition sample indices (for vertical markers) ─────────────
k_step1  = 101;
k_step2  = 301;
step_locs = [k_step1-1, k_step2-1];

%% ════════════════════════════════════════════════════════════════════════
%                     PERFORMANCE INDICES
% ════════════════════════════════════════════════════════════════════════

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
fprintf('           EXPERIMENT PERFORMANCE INDEX COMPARISON\n');
fprintf('           (Ns = %d samples,  %.1f min)\n', Ns, Ns*Ts/60);
fprintf('=================================================================\n');
fprintf('  Index          |   PI (Loop 1) |  PI (Loop 2) | MPC (Loop 1) | MPC (Loop 2)\n');
fprintf('-----------------+---------------+--------------+--------------+-------------\n');
fprintf('  SSE            |  %11.4f  |  %10.4f  |  %10.4f  | %10.4f\n', ...
    SSE_pi_1,  SSE_pi_2,  SSE_mpc_1,  SSE_mpc_2);
fprintf('  SSMV           |  %11.4f  |  %10.4f  |  %10.4f  | %10.4f\n', ...
    SSMV_pi_1, SSMV_pi_2, SSMV_mpc_1, SSMV_mpc_2);
fprintf('  SSdeltaMV      |  %11.4f  |  %10.4f  |  %10.4f  | %10.4f\n', ...
    SSdelMV_pi_1, SSdelMV_pi_2, SSdelMV_mpc_1, SSdelMV_mpc_2);
fprintf('=================================================================\n\n');

%% ════════════════════════════════════════════════════════════════════════
%                          COMPARISON PLOTS
% ════════════════════════════════════════════════════════════════════════

%% ── Figure 1: Yi(k) vs k and Ri(k) vs k ─────────────────────────────────
figure('Name','Exp Comparison – Output Profiles','NumberTitle','off','Position',[50 50 1000 700])

subplot(2,1,1)
plot(kT_cmp, Yk_pi(1,1:Ns),   'b-',  'LineWidth',1.8), hold on
plot(kT_cmp, Yk_mpc(1,1:Ns),  'r-',  'LineWidth',1.8)
plot(kT_cmp, Rk_f_pi(1,1:Ns), 'k--', 'LineWidth',1.2)
xline(step_locs(1),'k:','Step 1','LabelVerticalAlignment','bottom','FontSize',8)
xline(step_locs(2),'k:','Step 2','LabelVerticalAlignment','bottom','FontSize',8)
grid on
ylabel('T_1  (deg C)','FontSize',11)
title('Output Y_1(k) and Reference R_1(k) – Experiment','FontSize',12)
legend('Y_1^{PI}(k)','Y_1^{MPC}(k)','R_1(k)','Location','Best','FontSize',9)

subplot(2,1,2)
plot(kT_cmp, Yk_pi(2,1:Ns),   'b-',  'LineWidth',1.8), hold on
plot(kT_cmp, Yk_mpc(2,1:Ns),  'r-',  'LineWidth',1.8)
plot(kT_cmp, Rk_f_pi(2,1:Ns), 'k--', 'LineWidth',1.2)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('T_2  (deg C)','FontSize',11)
xlabel('Sample  k','FontSize',11)
title('Output Y_2(k) and Reference R_2(k)','FontSize',12)
legend('Y_2^{PI}(k)','Y_2^{MPC}(k)','R_2(k)','Location','Best','FontSize',9)

%% ── Figure 2: Ui(k) vs k ─────────────────────────────────────────────────
figure('Name','Exp Comparison – Heater Inputs','NumberTitle','off','Position',[100 50 1000 700])

subplot(2,1,1)
stairs(kT_cmp, Uk_pi(1,1:Ns)',  'b-', 'LineWidth',1.8), hold on
stairs(kT_cmp, Uk_mpc(1,1:Ns)', 'r-', 'LineWidth',1.8)
yline(5,  'k--', 'LineWidth',0.9)
yline(80, 'k--', 'LineWidth',0.9)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylim([0 105])
ylabel('U_1(k)  (% heater)','FontSize',11)
title('Heater Input U_1(k) – Experiment','FontSize',12)
legend('U_1^{PI}(k)','U_1^{MPC}(k)','Safety bounds','Location','Best','FontSize',9)

subplot(2,1,2)
stairs(kT_cmp, Uk_pi(2,1:Ns)',  'b-', 'LineWidth',1.8), hold on
stairs(kT_cmp, Uk_mpc(2,1:Ns)', 'r-', 'LineWidth',1.8)
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
figure('Name','Exp Comparison – Filtered Error','NumberTitle','off','Position',[150 50 1000 700])

subplot(2,1,1)
stairs(kT_cmp, ek_f_pi(1,1:Ns)',  'b-', 'LineWidth',1.8), hold on
stairs(kT_cmp, ek_f_mpc(1,1:Ns)', 'r-', 'LineWidth',1.8)
yline(0,'k--','LineWidth',0.8)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('e_{f,1}(k)  (deg C)','FontSize',11)
title('Filtered Error e_{f,1}(k) – Experiment','FontSize',12)
legend('PI','MPC','Location','Best','FontSize',9)

subplot(2,1,2)
stairs(kT_cmp, ek_f_pi(2,1:Ns)',  'b-', 'LineWidth',1.8), hold on
stairs(kT_cmp, ek_f_mpc(2,1:Ns)', 'r-', 'LineWidth',1.8)
yline(0,'k--','LineWidth',0.8)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('e_{f,2}(k)  (deg C)','FontSize',11)
xlabel('Sample  k','FontSize',11)
title('Filtered Error e_{f,2}(k)','FontSize',12)
legend('PI','MPC','Location','Best','FontSize',9)

%% ── Figure 4: xs(k) and us(k) vs k  (MPC only) ──────────────────────────
[n_st_bb, ~] = size(Xs_mpc);

figure('Name','Exp MPC – SS Targets xs(k) and us(k)','NumberTitle','off','Position',[200 50 1000 700])

subplot(2,1,1)
plot(kT_mpc(1:Ns), Xs_mpc(1,1:Ns), 'b-', 'LineWidth',1.8), hold on
if n_st_bb > 1
    plot(kT_mpc(1:Ns), Xs_mpc(2,1:Ns), 'r-', 'LineWidth',1.8)
end
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('x_{s,i}(k)  (deg C)','FontSize',11)
title('MPC SS Target States x_s(k)  [Experiment]','FontSize',12)
leg_xs = arrayfun(@(i) sprintf('x_{s,%d}^{MPC}(k)',i), 1:n_st_bb, 'UniformOutput',false);
legend(leg_xs{:},'Location','Best','FontSize',9)

subplot(2,1,2)
stairs(kT_mpc(1:Ns), Us_mpc(1,1:Ns)', 'b-', 'LineWidth',1.8), hold on
stairs(kT_mpc(1:Ns), Us_mpc(2,1:Ns)', 'r-', 'LineWidth',1.8)
xline(step_locs(1),'k:')
xline(step_locs(2),'k:')
grid on
ylabel('u_{s,i}(k)  (% heater)','FontSize',11)
xlabel('Sample  k','FontSize',11)
title('MPC SS Target Inputs u_s(k)','FontSize',12)
legend('u_{s,1}^{MPC}(k)','u_{s,2}^{MPC}(k)','Location','Best','FontSize',9)

%% ── Figure 5: Performance index bar chart ───────────────────────────────
figure('Name','Exp Performance Index Comparison','NumberTitle','off','Position',[250 50 1000 500])

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
sgtitle('Performance Indices – Experiment Comparison','FontSize',13)

%% ── Save ─────────────────────────────────────────────────────────────────
exp_compare.kT_cmp       = kT_cmp;
exp_compare.Ns            = Ns;
exp_compare.Yk_pi         = Yk_pi(:,1:Ns);
exp_compare.Uk_pi         = Uk_pi(:,1:Ns);
exp_compare.Rk_f_pi       = Rk_f_pi(:,1:Ns);
exp_compare.ek_f_pi       = ek_f_pi(:,1:Ns);
exp_compare.Yk_mpc        = Yk_mpc(:,1:Ns);
exp_compare.Uk_mpc        = Uk_mpc(:,1:Ns);
exp_compare.Rk_f_mpc      = Rk_f_mpc(:,1:Ns);
exp_compare.ek_f_mpc      = ek_f_mpc(:,1:Ns);
exp_compare.Xs_mpc        = Xs_mpc(:,1:Ns);
exp_compare.Us_mpc        = Us_mpc(:,1:Ns);
exp_compare.SSE_pi        = [SSE_pi_1;  SSE_pi_2 ];
exp_compare.SSE_mpc       = [SSE_mpc_1; SSE_mpc_2];
exp_compare.SSMV_pi       = [SSMV_pi_1;  SSMV_pi_2 ];
exp_compare.SSMV_mpc      = [SSMV_mpc_1; SSMV_mpc_2];
exp_compare.SSdelMV_pi    = [SSdelMV_pi_1;  SSdelMV_pi_2 ];
exp_compare.SSdelMV_mpc   = [SSdelMV_mpc_1; SSdelMV_mpc_2];

save TCL_Compare_ExpResults  exp_compare
fprintf('Results saved to  TCL_Compare_ExpResults.mat\n');
fprintf('All done.\n');
