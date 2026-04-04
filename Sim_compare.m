
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
     % change at k = 400
     if ( k < 250 )
        rk=[0;0] ;
     elseif (( k>= 250) && ( k < 450) )
        rk=[8;5] ;
     else
        rk=[0;0] ;
     end  

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

% Loading Experimental Data

if (linmod_type == 1)
    exp = load ('LQOC_BlackBox_Results.mat');
else
    exp = load ('LQOC_MechMod_Results.mat');
end

k_range=200:450;
N_use=length(k_range)-1;

% Taking only required experimental data

T1_exp = exp.t1s(k_range);
T2_exp = exp.t2s(k_range) ;
H1_exp = exp.h1s(k_range) ;  
H2_exp = exp.h2s(k_range) ;
R1_exp = exp.R1s(k_range) ;   
R2_exp = exp.R2s(k_range) ;

y1_exp = T1_exp - Ys(1) ;
y2_exp = T2_exp - Ys(2) ;
u1_exp = H1_exp - Us(1) ;
u2_exp = H2_exp - Us(2) ;
r1_exp = R1_exp - Ys(1) ;
r2_exp = R2_exp - Ys(2) ;

ef1_exp = exp.e1s(k_range) ;
ef2_exp = exp.e2s(k_range) ;

LQR_exp = exp.LQR;
xs1_exp = LQR_exp.xs_k(1, k_range) ;
xs2_exp = LQR_exp.xs_k(2, k_range) ;
us1_exp = LQR_exp.us_k(1, k_range) ;
us2_exp = LQR_exp.us_k(2, k_range) ;

% Taking only required Simulation data

y1_sim  = yk(1, k_range) ;
y2_sim  = yk(2, k_range) ;
u1_sim  = uk(1, k_range) ;
u2_sim  = uk(2, k_range) ;
r1_sim  = LQR.rk_f(1, k_range) ; 
r2_sim  = LQR.rk_f(2, k_range) ;
 
ef1_sim = LQR.ek_f(1, k_range) ;
ef2_sim = LQR.ek_f(2, k_range) ;
 
xs1_sim = LQR.xs_k(1, k_range) ;
xs2_sim = LQR.xs_k(2, k_range) ;
us1_sim = LQR.us_k(1, k_range) ;
us2_sim = LQR.us_k(2, k_range) ;

kT_plot_comp = k_range(1:N_use);


% Display simulation results

Set_Graphics_Parameters

Xk = TCL.Xs + xk ;   Uk = TCL.Us + uk ; 
Yk = yk + TCL.C_mat * TCL.Xs ;  Rk_f = LQR.rk_f + TCL.C_mat * TCL.Xs ;
u_limit = [LQR.u_L(1)*ones(N,1) LQR.u_H(1)*ones(N,1)]  ;


% Controlled output profiles 
figure
subplot(2,1,1)
plot(kT_plot_comp, y1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_comp, y1_exp(1:N_use), 'r--', 'LineWidth', 1.5)
plot(kT_plot_comp, r1_sim(1:N_use), 'k-.', 'LineWidth', 1)
grid on,
ylabel('y_1(k)')
title('Comparison: Output y_1 (Simulation vs Experiment)')
legend('Sim','Exp','Setpoint','Location', 'Best')
 
subplot(2,1,2)
plot(kT_plot_comp, y2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_comp, y2_exp(1:N_use), 'r--', 'LineWidth', 1.5)
plot(kT_plot_comp, r2_sim(1:N_use), 'k-.', 'LineWidth', 1)
grid on,
xlabel('Sampling Instant (k)'), ylabel('y_2(k)')
title('Comparison: Output y_2 (Simulation vs Experiment)')
legend('Sim','Exp','Setpoint','Location', 'Best')

% Manipulated Inputs u1 and u2 (perturbation)
figure
subplot(2,1,1)
stairs(kT_plot_comp, u1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
stairs(kT_plot_comp, u1_exp(1:N_use), 'r--', 'LineWidth', 1.5), grid on
ylabel('u_1(k)')
title('Comparison: Manipulated Input u_1 (Simulation vs Experiment)')
legend('Sim', 'Exp', 'Location', 'Best')
 
subplot(2,1,2)
stairs(kT_plot_comp, u2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
stairs(kT_plot_comp, u2_exp(1:N_use), 'r--', 'LineWidth', 1.5), grid on
xlabel('Sampling Instant (k)'), ylabel('u_2(k)')
title('Comparison: Manipulated Input u_2 (Simulation vs Experiment)')
legend('Sim', 'Exp', 'Location', 'Best')

% Target States xs(k)
figure
subplot(2,1,1)
plot(kT_plot_comp, xs1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_comp, xs1_exp(1:N_use), 'r--', 'LineWidth', 1.5), grid on
ylabel('x_{s,1}(k)'), 
title('Comparison: Target State x_{s,1}(k) (Simulation vs Experiment)')
legend('Sim', 'Exp', 'Location', 'Best')
 
subplot(2,1,2)
plot(kT_plot_comp, xs2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_comp, xs2_exp(1:N_use), 'r--', 'LineWidth', 1.5), grid on
xlabel('Sampling Instant (k)'), ylabel('x_{s,2}(k)')
title('Comparison: Target State x_{s,2}(k) (Simulation vs Experiment)')
legend('Sim', 'Exp', 'Location', 'Best')
 
% Target Inputs us(k)
figure
subplot(2,1,1)
stairs(kT_plot_comp, us1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
stairs(kT_plot_comp, us1_exp(1:N_use), 'r--', 'LineWidth', 1.5), grid on
ylabel('u_{s,1}(k)')
title('Comparison: Target Input u_{s,1}(k) (Simulation vs Experiment)')
legend('Sim', 'Exp', 'Location', 'Best')
 
subplot(2,1,2)
stairs(kT_plot_comp, us2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
stairs(kT_plot_comp, us2_exp(1:N_use), 'r--', 'LineWidth', 1.5), grid on
xlabel('Sampling Instant (k)'), ylabel('u_{s,2}(k)')
title('Comparison: Target Input u_{s,2}(k) (Simulation vs Experiment)')
legend('Sim', 'Exp', 'Location', 'Best')
 
% Filtered Innovation Sequence ef(k)
figure
subplot(2,1,1)
plot(kT_plot_comp, ef1_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_comp, ef1_exp(1:N_use), 'r--', 'LineWidth', 1.5), grid on
ylabel('e_{f,1}(k)')
title('Comparison: Filtered Innovation e_{f,1}(k) (Simulation vs Experiment)')
legend('Sim', 'Exp', 'Location', 'Best')
 
subplot(2,1,2)
plot(kT_plot_comp, ef2_sim(1:N_use), 'b-',  'LineWidth', 1.5), hold on
plot(kT_plot_comp, ef2_exp(1:N_use), 'r--', 'LineWidth', 1.5), grid on
xlabel('Sampling Instant (k)'), ylabel('e_{f,2}(k)')
title('Comparison: Filtered Innovation e_{f,2}(k) (Simulation vs Experiment)')
legend('Sim', 'Exp', 'Location', 'Best')
