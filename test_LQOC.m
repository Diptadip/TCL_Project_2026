close all; clear all; clc

% include tclab.m for initialization

tclab_N;

% Load LQOC Parameters 

N_samp = 450 ;

disp('Load Controller Files')

linmod_type = 2 ;  % Type = 1 for blackbox and 2 for mechanistic 
if (linmod_type == 1)
    load TCL_LQOC_BlackBox_Parameters
else
    load TCL_LQOC_MechMod_Parameters
end

n_ip = 2 ;  n_op = 2 ; 

LQR.xkhat = zeros(LQR.n_st,N_samp);  % LQR internal model states 
LQR.xs_k = zeros(LQR.n_st,N_samp);   % Target state 
LQR.us_k = zeros(n_ip,N_samp) ;        % Target input 

figure(1)
t1s = [];  t2s = [];
h1s = []; h2s = [];
e1s = [] ; e2s = [] ;
R1s = [] ; R2s = [] ;
Tk = zeros(n_op, 1) ; 

for k = 1:N_samp
    tic
    % read temperatures
    Tk(1) = T1C();
    Tk(2) = T2C();
    
    yk = Tk - Ys ;  % Find perturbation measurement  

    if ( k <= 50 )
         uk = [ 0 ; 0 ] ;   % Manual mode operation 
    else
        % 
        if ( k < 250 )
            setpt = zeros(n_op,1) ;
        elseif (( k>= 250) && ( k < 450) )
            setpt = [8 5]' ;
        else
             setpt = zeros(n_op,1) ;
        end
        
        % Implement MPC 

         [uk, LQR] = LQOC( LQR, setpt, yk, k )  ;
                
    end 
    
    Uk = uk + Us   ;
    Rk = LQR.rk_f(:,k) + Ys ; 
          
    h1(Uk(1));  % Send  manipulated inputs to heaters 
    h2(Uk(2));
        
     ek = LQR.ek_f(:,k) ;

    % plot heater and temperature data
    h1s = [h1s,Uk(1)];
    h2s = [h2s,Uk(2)];
    t1s = [t1s,Tk(1)];
    t2s = [t2s,Tk(2)];
    R1s = [ R1s, Rk(1) ] ;
    R2s = [ R2s, Rk(2) ] ;
    e1s = [e1s,ek(1)];
    e2s = [e2s,ek(2)];
    
    if ( rem(k,15)==0)  % Update plot after every 15 seconds 
        n = length(t1s);
        time = linspace(0,n+1,n);
        clf
        subplot(2,1,1)
        plot(time,t1s,'r-', time,R1s,'k.-','MarkerSize',10);
        hold on
        plot(time,t2s,'b-', time,R2s,'k.-','MarkerSize',10);
        ylabel('Temperature (degC)')
        legend('Temperature 1','Temperature 2','Location','NorthWest')
        subplot(2,1,2)
        plot(time,h1s,'r-', 'LineWidth',2);
        hold on
        plot(time,h2s,'b-','LineWidth',2);
        ylabel('Heater (0-5.5 V)')
        xlabel('Time (sec)')
        legend('Heater 1','Heater 2','Location','NorthWest')
        drawnow;
    end
    t = toc;
    [k LQR.samp_T-t] 
    pause(max(0.01,LQR.samp_T-t))
end

disp('Turn off heaters')
h1(0);
h2(0);

disp('LQOC Experimental Test Complete')

% Plot experimental results 

tf = N_samp * LQR.samp_T ;

figure(2), subplot(2,1,1) 
plot(time,t1s,'r-', time,R1s,'k.-','MarkerSize',10), grid
axis([0 tf 30 70])
ylabel('Temperature (degC)')
legend('Temperature 1','Setpoint 1')
subplot(2,1,2)
plot(time,t2s,'b-', time,R2s,'k.-','MarkerSize',10), grid
axis([0 tf 30 70]) 
legend('Temperature 2','Setpoint 2')
ylabel('Temperature (degC)')
xlabel('Time (sec)')

figure(3), subplot(2,1,1),
stairs(time,h1s,'r-', 'LineWidth',2), grid
axis([0 tf 0 100]) 
ylabel('Heater 1 (%)')
subplot(2,1,2), stairs(time,h2s,'b-','LineWidth',2), grid
axis([0 tf 0 100]) 
ylabel('Heater 2 (%)')
xlabel('Time (sec)')

figure(4), subplot(2,1,1) 
plot(time,e1s,'r-','MarkerSize',10), grid
axis([0 tf -5 5]) 
title('Filtered Innovation Sequence') 
ylabel('ek_f_,_1(k)')
subplot(2,1,2) 
plot(time,e2s,'b-','MarkerSize',10), grid
axis([0 tf -5 5]) 
ylabel('ek_f_,_2(k)')
xlabel('Time (sec)')
      
% Save experimental data 

if (linmod_type == 1)
   save LQOC_BlackBox_Results time h1s h2s t1s t2s R1s R2s e1s e2s LQR
else
    save LQOC_MechMod_Results time h1s h2s t1s t2s R1s R2s e1s e2s LQR
end



