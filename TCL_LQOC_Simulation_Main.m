
clear all
close all 

global TCL 

load TCL_MechModel_Parameters.mat
load TCL_BlackBox_OE_SS.mat

% Initialization for Process Simulation 
N = 400 ;   % Set Simulation time 

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
