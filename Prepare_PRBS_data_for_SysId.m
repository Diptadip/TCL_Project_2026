% This file prepares perturbation data sets for model identification and
% validation. It also estimates operating steady state (Ys,Us) and
% measurement noise matrix R from the experimental data 

clear all
close all
clc

load Data\PRBS_Data_Group9_23march26.mat

% SetGraphics

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

% Estimate the output steady state using 150 samples before PRBS starts 
N0 = 850/samp_T ;              % Output steady state 
N1 = 1000/samp_T ; 

figure, subplot(211), plot( time(N0:N1), Yk_data(N0:N1,1) ), grid 
ylabel('T_1(k)'), title('Raw Output Data' )
subplot(212), plot( time(N0:N1), Yk_data(N0:N1,2) ), grid 
ylabel('T_2(k)'), xlabel('Time (sec.)')

Ys = mean( Yk_data(N0:N1,:) ) 
R_mat = diag(var( Yk_data(N0:N1,:) ))    % Measurement noise covariance 

% Generate perturbation data for system identification 

time0 = time ;  % Save the initial time vector 
clear time 

N1 = 3600 / samp_T ;    % Redefine N1 to end of identification data set

yk_id = Yk_data(N0:N1,:) - Ys ;
uk_id = Uk_data(N0:N1,:) - Us' ;
time = time0(N0:N1) ;
id_data = iddata( yk_id, uk_id, samp_T ) ;    % Create identification data set 

% Plot identification data 

figure, subplot(211), plot( time, yk_id(:,1) ), grid 
ylabel('y_1(k)'), title('Perturbation Output Data' )
subplot(212), plot( time, yk_id(:,2) ), grid 
ylabel('y_2(k)'), xlabel('Time (sec.)')

figure,subplot(211), stairs( time, uk_id(:,1), 'LineWidth',2 ), grid 
title('Perturbation Input Data' ), ylabel('u_1(k)') 
subplot(212), stairs( time, uk_id(:,2), 'LineWidth',2 ), grid 
xlabel('Time (sec.)'), ylabel('u_2(k)') 

% Create validation data set 

N1 = 3000 / samp_T ;  % Reset N1 to beginning of validation data 
yk_val = Yk_data(N1:end,:) - Ys ;
uk_val = Uk_data(N1:end,:) - Us' ;

clear time
time = time0(N1:end) ;
val_data = iddata( yk_val, uk_val, samp_T ) ;

% Plot validation data 

figure, subplot(211), plot( time, yk_val(:,1) ), grid 
ylabel('y_1(k)'), title('Perturbation Output Data' )
subplot(212), plot( time, yk_val(:,2) ), grid 
ylabel('y_2(k)'), xlabel('Time (sec.)')

figure,subplot(211), stairs( time, uk_val(:,1), 'LineWidth',2 ), grid 
title('Perturbation Input Data' ), ylabel('u_1(k)') 
subplot(212), stairs( time, uk_val(:,2), 'LineWidth',2 ), grid 
xlabel('Time (sec.)'), ylabel('u_2(k)') 

save SysId_Data_Sets id_data val_data Ys Us R_mat
 





