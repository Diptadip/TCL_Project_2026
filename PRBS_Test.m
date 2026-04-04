% This program injects PRBS inputs in TC Lab system 
% System is brought to a steady state in the beginning and then
% PRBS perturbations are introduced 

close all; clear all; clc

% include tclab.m for initialization
tclab_N;

n_ip = 2 ;  % Number of inputs 

disp('Turn off heaters')
h1(0);
h2(0);

disp('Test Heater 1')
disp('LED Indicates Temperature')

figure(1)

tfinal = 4200 ;  % Time in seconds
% tfinal = 200 ;  % Used for program testing
samp_T = 4 ; % Ts = 4 sec 
N_samp = tfinal / samp_T ;

% Define arrays for saving data 
T1_s = zeros(N_samp,1); 
T2_s = zeros(N_samp,1) ;
H1_s = zeros(N_samp,1);
H2_s = zeros(N_samp,1);

% PRBS parameters 
N_start= 1000 / samp_T ;          % PRBS start time 
% N_start = 10 ;    % Used for program testing
switch_period = 48 / samp_T ;  % Switching period 
prbs_ampl = [ 10 ; 12 ]                          % PRBS amplitude 
Hs = [ 30 ; 40 ] ;                         % Initial steady state input 

time = zeros(N_samp,1) ;
display_time = 16 / samp_T ;

% initial heater values

% Perform PRBS test 
for k = 1:N_samp
    tic;
    time(k) = (k-1) * samp_T ; 
    
    % Get temperature measurements from TC Lab  
    Tk(1) = T1C();
    Tk(2) = T2C();
    
    % Input signal generation 
    if ( k < N_start)
        delta_H = zeros(n_ip,1) ;   % Move system to steady state at Hk = Hs 
    else
        if (rem(k, switch_period )==0)  % Check PRBS switching instant 
             delta_H = Gen_PRBS_Signal(  n_ip, prbs_ampl ) ; %  PRBS generator
        end
    end
    
    % Send manipulated inputs to TC lab 
    Hk = Hs + delta_H ;
    h2(Hk(2));   % Inject input 2 into the system
    h1(Hk(1));   % Inject input 1 into the system
    
    % Save input and output data 
     H1_s(k) = Hk(1);      H2_s(k) = Hk(2);
     T1_s(k) = Tk(1);      T2_s(k) = Tk(2);
    
    % plot heater and temperature data
    if (rem(k,display_time )==0)  % Generate random numbers
        clf
        subplot(2,1,1)
        plot(time(1:k),T1_s(1:k),'r.','MarkerSize',8);
        hold on
        plot(time(1:k),T2_s(1:k),'b.','MarkerSize',8);
        ylabel('Temperature (degC)')
        legend('Temperature 1','Temperature 2','Location','NorthWest')
        subplot(2,1,2)
        stairs(time(1:k),H1_s(1:k),'r-','LineWidth',2);
        hold on
        stairs(time(1:k),H2_s(1:k),'b--','LineWidth',2);
        ylabel('Heater (0-100%)');xlabel('Time (sec)')
        legend('Heater 1','Heater 2','Location','NorthWest')
        drawnow;
    end
    t = toc;
    [k samp_T-t] 
    pause(max(0.01, samp_T -t))   % Wait till sampling interval is over 
end

disp('Turn off heaters')
h1(0);
h2(0);

% Plot results 

figure(2), subplot(2,1,1) 
plot(time,T1_s,'r-'), grid
axis([0 N_samp*samp_T 30 70])
ylabel('Temperature (degC)')
subplot(2,1,2)
plot(time,T2_s,'b-'), grid
axis([0 N_samp*samp_T 30 70]) 
ylabel('Temperature (degC)')
xlabel('Time (sec)')

figure(3), subplot(2,1,1),
stairs(time,H1_s,'r-', 'LineWidth',2), grid
axis([0 N_samp*samp_T 30 70]) 
ylabel('Heater 1 (%)')
subplot(2,1,2), stairs(time,H2_s,'b-','LineWidth',2), grid
axis([0 N_samp*samp_T 30 70]) 
ylabel('Heater 2 (%)')
xlabel('Time (sec)')

% Save expt data 
Us = Hs ; 
Uk_data = [ H1_s H2_s ] ;
Yk_data = [ T1_s T2_s ] ; 

save PRBS_Data_GroupX_Date time Yk_data Uk_data Us  samp_T

disp('Heater Test Complete')

   