
clear all
% close all 

global TCL 

load TCL_MechModel_Parameters.mat
load TCL_BlackBox_OE_SS.mat

% Initialization for Process Simulation 
N = 450 ;   % Set Simulation time 

[n_op, nx ] = size(TCL.C_mat)  ;
[n_st, n_ip] = size(TCL.gama) ;

Ys = TCL.Ys ;  Us = TCL.Us ;    % These are needed in experimental implementation

% Initialize internal model for LQOC

linmod_type = 2 ;  % Type = 1 for blackbox and 2 for mechanistic 

if (linmod_type == 1)   % Blackbox model 
    LQR.phy=idmod.phy;
    LQR.gam=idmod.gama;
    LQR.C_mat=idmod.C_mat;
    LQR.L_mat=idmod.L_mat;
else                                       % Linearized mechanistis model 
    LQR.phy=TCL.phy;
    LQR.gam= TCL.gama;
    LQR.C_mat= TCL.C_mat;  % Temperature T(1) and T(2) are measured in TCL
    LQR.L_mat= zeros(n_st,n_op);
end
LQR.samp_T = TCL.Samp_T ;

[LQR.n_st, LQR.n_ip] = size(LQR.gam) ;
LQR.n_op = n_op ;

LQR.u_H = [100 100]' - TCL.Us ;  % Physical Input bounds %-heating 
LQR.u_L = [0 0]' - TCL.Us ;  

% Using these tuning parameters and the model obtained from system identification, controller gain is obtained using dlqr command in MATLAB
if (linmod_type == 1)
    beta(1)=1; beta(2)=1; alpha(1)=50; alpha(2)=50;
else
    beta(1)=1; beta(2)=1; alpha(1)=10; alpha(2)=10;
end

% LQOC tuning matrices and controller design 
LQR.Wu = diag(beta) ;
LQR.Wy =diag(alpha) ;

Wx_y = LQR.C_mat' * LQR.Wy * LQR.C_mat ;
[LQR.G_inf,S_inf, Eig_CL_con] = dlqr(LQR.phy, LQR.gam, Wx_y, LQR.Wu, zeros(LQR.n_st,LQR.n_ip)) ;
disp('Controller gain is')
LQR.G_inf
disp('Maximum Absolute Closed loop Eigen value')
rho_con = max(abs(Eig_CL_con)) 

% Gains used in the computation of target states and target inputs.

LQR.Gain_ss = LQR.C_mat * inv(eye(LQR.n_st)-LQR.phy)*LQR.gam ; 
LQR.Gain_ss_inv =inv(LQR.Gain_ss) ;
LQR.inv_mat = inv(eye(LQR.n_st) - LQR.phy) ;
LQR.Gain_e = eye(LQR.n_op) + LQR.C_mat * inv(eye(LQR.n_st)-LQR.phy)* LQR.L_mat ;

alfa_e =0.95;  % Set this value to 0 to see effect of noise on u(k)
LQR.phy_e = alfa_e * eye(n_op) ;  % Innovation filter
LQR.ek = zeros(n_op,N) ;
LQR.ek_f = zeros(n_op,N) ;

% Task 5: Change setpoint filter parameter between 0.9 and 0.99. Setting it equal to zero will eliminate innovation filter. 
beta_r =0.95;  % Set this value to 0 see effect of step change in rk on u(k)
LQR.phy_r = beta_r * eye(n_op) ; % 
LQR.rk_f = zeros(n_op,N);     % Setpoint filter parameters

if (linmod_type == 1)
    save TCL_LQOC_BlackBox_Parameters LQR Ys Us
else
    save TCL_LQOC_MechMod_Parameters LQR Ys Us
end

LQR.xkhat = zeros(LQR.n_st,N);  % LQR internal model states 
LQR.xs_k = zeros(LQR.n_st,N);   % Target state 
LQR.us_k = zeros(n_ip,N) ;        % Target input 

xk = zeros(n_st,N);   % Simulating linear mechanistic plant 
uk = zeros(n_ip,N) ;
yk = zeros(n_op,N) ;
yk(:,1) = TCL.C_mat * xk(:,1) ;

Noise_ON = 1 ; 
vk = mvnrnd(zeros(n_op,1), TCL.R, N) ;
vk = Noise_ON * vk' ;

dk = zeros(1,N) ;

rk = zeros(n_op,1) ;

% Select Plant_Sim as Linear for simulating Q Tank process using linear discrete model.

Plant_Sim =1  ;

for k = 1:N-1
    k
    kT(k) = (k-1);

     % Setpoint change introduced simultaneously with the disturbance
     % change at k = 250
     if (k >= 50),   rk = [ 8 5 ]';  end

      [uk_LQR, LQR] = LQOC( LQR, rk, yk(:,k), k ) ;
       uk(:,k) = uk_LQR ;

     % Quadruple tank Plant Simulation from k to k+1

    if (Plant_Sim )
        xk(:,k+1) = TCL.phy * xk(:,k) + TCL.gama * uk(:,k) ;
    else
        TCL.Uk = TCL.Us + uk(:,k) ;  
        % Nonlinear plant simulation
        [T,Xt] = ode45( 'TCL_Dynamics' , [ 0 TCL.Samp_T ], TCL.Xs+xk(:,k)) ;  
        xk(:,k+1) = Xt( end,: )' - TCL.Xs ;
    end
    if (Noise_ON)
        yk(:,k+1) = TCL.C_mat * xk(:,k+1) + vk(:,k+1);  %Measurement with noise
    else
        yk(:,k+1) = TCL.C_mat * xk(:,k+1) ;  %Measurement with noise
    end

 end
uk(:,N) = uk(:,N-1) ;
LQR.ek(:,N) = LQR.ek(:,N-1) ;
LQR.xs_k(:,N) = LQR.xs_k(:,N-1) ;
LQR.us_k(:,N) = LQR.us_k(:,N-1) ; 
kT(N) = N-1 ;

% exp = load('LQOC_MechMod_Results.mat');
% idx=200:450;
% time_idx=exp.time(idx);

% Display simulation results

Set_Graphics_Parameters

Xk = TCL.Xs + xk ;   Uk = TCL.Us + uk ; 
Yk = yk + TCL.C_mat * TCL.Xs ;  Rk_f = LQR.rk_f + TCL.C_mat * TCL.Xs ;
u_limit = [LQR.u_L(1)*ones(N,1) LQR.u_H(1)*ones(N,1)]  ;

%Controlled output profiles 
figure, subplot(211), 
plot(kT, yk(1,:), 'b',  kT, LQR.rk_f(1,:), '-.', kT, yk(2,:), 'k', kT, LQR.rk_f(2,:), '-.'), grid
title('Controlled Output and Setpoint Profiles (Deviation)'), ylabel('y_i (k) (volts)')
subplot(212), plot(kT, Yk(1,:), 'b',  kT, Rk_f(1,:), '-.', kT, Yk(2,:), 'k', kT, Rk_f(2,:), '-.'), grid
xlabel('Sampling Instant (k)'), ylabel('Y_i (k) (volts)')
legend('h_1 (k)', 'r_1 (k)', 'h_2 (k)', 'r_2 (k)', 'Orientation','horizontal','Location','Best')
title('Controlled Output and Setpoint Profiles')

% Manipulated input profiles
figure, subplot(211), stairs(kT, LQR.us_k(1,:)', '-.b','LineWidth',1), grid
hold on
stairs(kT, uk(1,:)', 'b','LineWidth',1)
stairs(kT, LQR.us_k(2,:)', '-.k','LineWidth',1)
stairs(kT, uk(2,:)', 'k','LineWidth',1)
stairs(u_limit, '-.r', 'LineWidth',1)
hold off
ylabel('u_i (k)'), title('Target Input Profiles and Deviation Input profiles')
legend( 'us_1', 'u_1', 'us_2', 'u_2', 'Orientation','horizontal','Location','Best')
subplot(212), stairs(kT, Uk(1,:)', 'b','LineWidth',1),  grid
hold on, stairs(kT, Uk(2,:)', 'k','LineWidth',1),
stairs(u_limit+TCL.Us', '-.r', 'LineWidth',1), hold off
ylabel('u_i (k)')
title('Manipulated Input Profiles')
legend( 'U_1 (k)', 'U_2 (k)','Orientation','horizontal','Location','Best')

% Disturbance input and Model Plant Mismatch profiles 
figure, stairs( kT, LQR.ek_f','LineWidth',1 ), grid
xlabel('Sampling Instant (k)'), ylabel('e_f_,_i (k)')
title('Model plant mismatch signal')

%System State profiles and LQR internal model profiles 

figure, plot(kT, Xk(1,:), kT, Xk(2,:)), grid
xlabel('Sampling Instant (k)'), ylabel('x_i (k)')
legend('T_1 (k)', 'T_2 (k)', 'Orientation','horizontal','Location','Best')
title('System State Profiles')

figure, plot(kT, LQR.xs_k, '-.r', kT, LQR.xkhat, 'k'), grid
title('RED: Target States   BLACK: LQR Model State Estimates')
xlabel('Sampling Instant (k)'), ylabel('xs_i(k) and x_i(k)')

k_range_sim = 200 : min(450, N) ;
kT_comp     = k_range_sim - 1 ;

y1_sim  = yk(1, k_range_sim) ;
y2_sim  = yk(2, k_range_sim) ;
u1_sim  = uk(1, k_range_sim) ;
u2_sim  = uk(2, k_range_sim) ;
r1_sim  = LQR.rk_f(1, k_range_sim) ;  % filtered setpoint (perturbation)
r2_sim  = LQR.rk_f(2, k_range_sim) ;
 
ef1_sim = LQR.ek_f(1, k_range_sim) ;
ef2_sim = LQR.ek_f(2, k_range_sim) ;
 
xs1_sim = LQR.xs_k(1, k_range_sim) ;
xs2_sim = LQR.xs_k(2, k_range_sim) ;
us1_sim = LQR.us_k(1, k_range_sim) ;
us2_sim = LQR.us_k(2, k_range_sim) ;

disp('Loading experimental data for comparison...')
if (linmod_type == 1)
    load LQOC_BlackBox_Results
else
    load LQOC_MechMod_Results
end

k_range_exp = 200:450 ;          % sample indices in experimental data
N_comp      = length(k_range_exp);  % = 251 points

% Absolute experimental temperatures and heater values
T1_exp = t1s(k_range_exp) ;    % (degC)
T2_exp = t2s(k_range_exp) ;
H1_exp = h1s(k_range_exp) ;    % (% heater, absolute)
H2_exp = h2s(k_range_exp) ;
R1_exp = R1s(k_range_exp) ;    % filtered setpoint (absolute)
R2_exp = R2s(k_range_exp) ;
 
% Convert to perturbation variables (subtract operating point)
y1_exp = T1_exp - Ys(1) ;
y2_exp = T2_exp - Ys(2) ;
u1_exp = H1_exp - Us(1) ;
u2_exp = H2_exp - Us(2) ;
r1_exp = R1_exp - Ys(1) ;
r2_exp = R2_exp - Ys(2) ;

% Filtered innovation sequence (already perturbation-form in test_LQOC.m)
ef1_exp = e1s(k_range_exp) ;
ef2_exp = e2s(k_range_exp) ;



% Retrieve xs_k and us_k from the experimental LQR structure
% (saved as LQR inside the experimental .mat file; rename to avoid overwriting sim LQR)
LQR_exp = LQR ;   % rename to avoid conflict
xs1_exp = LQR_exp.xs_k(1, k_range_exp) ;
xs2_exp = LQR_exp.xs_k(2, k_range_exp) ;
us1_exp = LQR_exp.us_k(1, k_range_exp) ;
us2_exp = LQR_exp.us_k(2, k_range_exp) ;

 



% Align lengths (in case N < 450, the simulation window is shorter)
N_use = min(length(k_range_sim), N_comp) ;
kT_plot_sim = kT_comp(1:N_use) ;
kT_plot_exp = k_range_exp(1:N_use) - 1 ;   % same x-axis

fprintf('\n========== LQOC Tuning Parameters ==========\n')
if (linmod_type == 1)
    fprintf('Model type       : Black-box (OE State-Space)\n')
else
    fprintf('Model type       : Linearized Mechanistic\n')
end
fprintf('alpha (Wy diag)  : [%.2f, %.2f]  (output weight)\n', alpha(1), alpha(2))
fprintf('beta  (Wu diag)  : [%.2f, %.2f]  (input weight)\n',  beta(1),  beta(2))
fprintf('alfa_e           : %.2f           (innovation filter pole)\n', alfa_e)
fprintf('beta_r           : %.2f           (setpoint filter pole)\n',   beta_r)
fprintf('Closed-loop rho  : %.4f\n', rho_con)
fprintf('=============================================\n\n')

% -----------------------------------------------------------------------
% COMPARISON FIGURE 1: Controlled Outputs Y1 and Y2 (perturbation)
% -----------------------------------------------------------------------
figure
subplot(2,1,1)
plot(kT_plot_sim, y1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_exp, y1_exp(1:N_use), 'r--', 'LineWidth', 1.5)
plot(kT_plot_sim, r1_sim(1:N_use), 'k-.', 'LineWidth', 1)
hold off, grid
ylabel('y_1(k)  [dev, degC]')
title('Comparison: Output y_1 (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Setpoint (Sim)', 'Location', 'Best')
 
subplot(2,1,2)
plot(kT_plot_sim, y2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_exp, y2_exp(1:N_use), 'r--', 'LineWidth', 1.5)
plot(kT_plot_sim, r2_sim(1:N_use), 'k-.', 'LineWidth', 1)
hold off, grid
xlabel('Sampling Instant (k)'), ylabel('y_2(k)  [dev, degC]')
title('Comparison: Output y_2 (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Setpoint (Sim)', 'Location', 'Best')
 
% -----------------------------------------------------------------------
% COMPARISON FIGURE 2: Manipulated Inputs U1 and U2 (perturbation)
% -----------------------------------------------------------------------
figure
subplot(2,1,1)
stairs(kT_plot_sim, u1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
stairs(kT_plot_exp, u1_exp(1:N_use), 'r--', 'LineWidth', 1.5)
hold off, grid
ylabel('u_1(k)  [dev, %]')
title('Comparison: Manipulated Input u_1 (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Location', 'Best')
 
subplot(2,1,2)
stairs(kT_plot_sim, u2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
stairs(kT_plot_exp, u2_exp(1:N_use), 'r--', 'LineWidth', 1.5)
hold off, grid
xlabel('Sampling Instant (k)'), ylabel('u_2(k)  [dev, %]')
title('Comparison: Manipulated Input u_2 (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Location', 'Best')
 
% -----------------------------------------------------------------------
% COMPARISON FIGURE 3: Target States xs(k)
% -----------------------------------------------------------------------
figure
subplot(2,1,1)
plot(kT_plot_sim, xs1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_exp, xs1_exp(1:N_use), 'r--', 'LineWidth', 1.5)
hold off, grid
ylabel('x_{s,1}(k)'), 
title('Comparison: Target State x_{s,1}(k) (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Location', 'Best')
 
subplot(2,1,2)
plot(kT_plot_sim, xs2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_exp, xs2_exp(1:N_use), 'r--', 'LineWidth', 1.5)
hold off, grid
xlabel('Sampling Instant (k)'), ylabel('x_{s,2}(k)')
title('Comparison: Target State x_{s,2}(k) (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Location', 'Best')
 
% -----------------------------------------------------------------------
% COMPARISON FIGURE 4: Target Inputs us(k)
% -----------------------------------------------------------------------
figure
subplot(2,1,1)
stairs(kT_plot_sim, us1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
stairs(kT_plot_exp, us1_exp(1:N_use), 'r--', 'LineWidth', 1.5)
hold off, grid
ylabel('u_{s,1}(k)  [dev, %]')
title('Comparison: Target Input u_{s,1}(k) (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Location', 'Best')
 
subplot(2,1,2)
stairs(kT_plot_sim, us2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
stairs(kT_plot_exp, us2_exp(1:N_use), 'r--', 'LineWidth', 1.5)
hold off, grid
xlabel('Sampling Instant (k)'), ylabel('u_{s,2}(k)  [dev, %]')
title('Comparison: Target Input u_{s,2}(k) (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Location', 'Best')
 
% -----------------------------------------------------------------------
% COMPARISON FIGURE 5: Filtered Innovation Sequence ef(k)
% -----------------------------------------------------------------------
figure
subplot(2,1,1)
plot(kT_plot_sim, ef1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_exp, ef1_exp(1:N_use), 'r--', 'LineWidth', 1.5)
hold off, grid
ylabel('e_{f,1}(k)')
title('Comparison: Filtered Innovation e_{f,1}(k) (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Location', 'Best')
 
subplot(2,1,2)
plot(kT_plot_sim, ef2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_exp, ef2_exp(1:N_use), 'r--', 'LineWidth', 1.5)
hold off, grid
xlabel('Sampling Instant (k)'), ylabel('e_{f,2}(k)')
title('Comparison: Filtered Innovation e_{f,2}(k) (Simulation vs Experiment) — Samples 200 to 450')
legend('Simulation', 'Experiment', 'Location', 'Best')
 
disp('Simulation vs Experiment comparison plots generated.')
