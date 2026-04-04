clc
clear all
close all

SetGraphics 

load SysId_Data_Sets.mat
load TCL_BlackBox_OE_SS    % Load identified black-box model 
load TCL_MechModel_Parameters  % Load TCL mechanistic model 

% Find simulated outputs using the identified black-box model 
uk = id_data.InputData' ;
yk = id_data.OutputData' ;

N = length(uk) ;
n_st = length( idmod.phy ) ;
n_op = 2 ; 

xk = zeros(2,N) ;   % Matrix for simulation  of linearized mechanistic state space model 
yk_hat = zeros(n_op,N) ;

zk = zeros(n_st,N) ;  % matrix for simulation of black-box state space model 
yk_id = zeros(n_op,N) ;

% Open loop Simulation 
for k = 1 : N-1
    % Simulate linear mechanistic model 
    xk(:,k+1) = TCL.phy * xk(:,k) + TCL.gama * uk(:,k) ; 
    yk_hat(:,k+1) = TCL.C_mat * xk(:,k+1) ;

    % Simulate linear black-box model 
    zk(:,k+1) = idmod.phy * zk(:,k) + idmod.gama * uk(:,k) ; 
    yk_id(:,k+1) = idmod.C_mat * zk(:,k+1) ;
end

figure, plot( yk_id(1,:), 'k')
hold on
plot( yk(1,:),'xr')
plot( yk_hat(1,:),'b')
xlabel('Sampling Instant'), ylabel('y_1(k)')
legend('Black-Box', 'Data', 'Mechanistic')
title('Identification Data Set')
hold off

figure, plot( yk_id(2,:), 'k')
hold on 
plot( yk(2,:),'xr')
plot( yk_hat(2,:),'b')
xlabel('Sampling Instant'), ylabel('y_2(k)')
legend('Black-Box', 'Data', 'Mechanistic')
title('Identification Data Set')
hold off

clear yk uk xk zk yk_hat yk_id 

uk = val_data.InputData' ;
yk = val_data.OutputData' ;

[n_ip, N] = size(uk) ;

xk = zeros(2,N) ;   % Matrix for simulation  of linearized mechanistic state space model 
yk_hat = zeros(n_op,N) ;

zk = zeros(n_st,N) ;  % matrix for simulation of black-box state space model 
yk_id = zeros(n_op,N) ;

% Open loop Simulation 
for k = 1 : N-1
    % Simulate linear mechanistic model 
    xk(:,k+1) = TCL.phy * xk(:,k) + TCL.gama * uk(:,k) ; 
    yk_hat(:,k+1) = TCL.C_mat * xk(:,k+1) ;

    % Simulate linear black-box model 
    zk(:,k+1) = idmod.phy * zk(:,k) + idmod.gama * uk(:,k) ; 
    yk_id(:,k+1) = idmod.C_mat * zk(:,k+1) ;
end

figure, plot( yk_id(1,:), 'k')
hold on
plot( yk(1,:),'xr')
plot( yk_hat(1,:),'b')
xlabel('Sampling Instant'), ylabel('y_1(k)')
legend('Black-Box', 'Data', 'Mechanistic')
title('Validation Data Set')
hold off

figure, plot( yk_id(2,:), 'k')
hold on 
plot( yk(2,:),'xr')
plot( yk_hat(2,:),'b')
xlabel('Sampling Instant'), ylabel('y_2(k)')
legend('Black-Box', 'Data', 'Mechanistic')
title('Validation Data Set')
hold off

% Compare model step responses and frequency responses 

id_mod = ss( idmod.phy, idmod.gama, idmod.C_mat, idmod.D_mat, idmod.samp_T) 
lin_mod = ss( TCL.phy, TCL.gama, TCL.C_mat, TCL.D_mat, TCL.Samp_T) 

ltiview( id_mod, lin_mod )