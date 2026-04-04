% This program identifies 3 heat transfer related parameters of the mechanistic model 
% for Temperature Control Lab system using the raw PRBD test experimental
% data set and the mechanistic model available in the literature. It also estimates operating steady 
% state (Ys,Us) and measurement noise matrix R from the experimental data 

clear all
close all
clc

load Data\PRBS_Data_Group9_23march26.mat

global TCL Yk_data  Uk_data Xk_data 

[nr, nc] = size(Yk_data) ;
Xk_data = zeros(nr,nc) ;

% Fixed parameters

TCL.eps = 0.9;                  % emissivity
TCL.sigma = 5.67*10^-8;  % Stefan Boltzman constant
TCL.m = 0.004;                % mass
TCL.As = 2*10^-4;            % surface area between heaters
TCL.A = 1*10^-3;              % surface area not between heaters
TCL.Cp = 500;                   % heat capacity
TCL.Tinf = 301;                % atmospheric temperature  
TCL.Samp_T = 4;             % Sampling interval is 4 sec

% The first step is to view the experimental data 

figure, subplot(211), plot( time, Yk_data(:,1) ), grid 
ylabel('T_1(k)'), title('Raw Output Data' )
subplot(212), plot( time, Yk_data(:,2) ), grid 
ylabel('T_2(k)'), xlabel('Time (sec.)')

figure,subplot(211), stairs( time, Uk_data(:,1), 'LineWidth',2 ), grid 
title('Raw Heater Input Data' ), ylabel('H_1(k)') 
subplot(212), stairs( time, Uk_data(:,2), 'LineWidth',2 ), grid 
xlabel('Time (sec.)'), ylabel('H_2(k)') 

% Find the data segment that can be used for estimating the steady state operating point (Ys, Us) 

% Estimate the output steady state and measurement noise 

N0 = 850/samp_T ;              % Output steady state 
N1 = 1000/samp_T ; 

figure, subplot(211), plot( time(N0:N1), Yk_data(N0:N1,1) ), grid 
ylabel('T_1(k)'), title('Raw Output Data' )
subplot(212), plot( time(N0:N1), Yk_data(N0:N1,2) ), grid 
ylabel('T_2(k)'), xlabel('Time (sec.)')

TCL.C_mat = eye(2)  ;    % TCL Measurement model 
TCL.D_mat = zeros(2) ; 

Ys = mean( Yk_data(N0:N1,:) ) ;
Us = mean( Uk_data(N0:N1,:) ) ;

TCL.R = diag(var( Yk_data(N0:N1,:) )) ;   
TCL.Xs = Ys' ;
TCL.Ys = TCL.C_mat * TCL.Xs ;
TCL.Us = Us' ; 

% Weighted least squares optimization for model parameter estimation 

% Transpose taken for the objective function calculations 
Yk_data  = Yk_data'+273; % Temperature converted to Kelvins from deg C    
Uk_data = Uk_data' ; 

% options = optimset('Display', 'iter', 'TolFun', 1e-4) ;
oldopts = optimset ;
tolerance = 1e-10 ;
options = optimset(oldopts,'MaxIter', 1e6, 'Display','iter', 'TolFun',tolerance, 'TolX', tolerance) ;

parameters0=[8, 0.005, 0.003]; % initial guess of model parameters
parameters = fmincon(@obj_paraEstm_TCL, parameters0, [], [], [], [], [0,0,0], [100,1,1],[],options) ; 

%  Transpose taken again for the purpose of column-wise plotting 
Yk_data  =Yk_data';   
Xk_data = Xk_data' ; 

figure, plot( time, Yk_data(:,1), time, Xk_data(:,1) ), grid 
ylabel('T_1(k)'), xlabel('Time (sec.)'),
title('Comparison of measured and estimated output 1' )

figure, plot( time, Yk_data(:,2), time, Xk_data(:,2) ), grid 
ylabel('T_2(k)'), xlabel('Time (sec.)')
title('Comparison of measured and estimated output 2' )

TCL.U = parameters(1);          % Overall Heat transfer coefficient
TCL.alpha(1) = parameters(2);  % Heater 1 efficiency parameter 
TCL.alpha(2) = parameters(3);  % Heater 2 efficiency parameter 

save TCL_MechModel_Parameters TCL


