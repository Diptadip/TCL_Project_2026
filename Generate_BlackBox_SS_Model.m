
%%%% Black-Box Output Error STATE SPACE Model IDENTIFICATION  %%%%%%
% using System Identification Toolbox of MATLAB 

clc
clear all
close all

load SysId_Data_Sets.mat

n_ip = 2 ; n_op = 2 ;  
n_st = input('Specify State Dimension (bet. 2 and 5): ')
BlackBox_model = pem(id_data, n_st, 'DisturbanceModel', 'none')

% Find simulated outputs using the identified black-box model 
uk = id_data.InputData ;
yk = id_data.OutputData ;

yk_id = idsim(BlackBox_model, uk);

% Save model matrices for later use in controller synthesis 
idmod.phy= BlackBox_model.A
idmod.gama = BlackBox_model.B 
idmod.C_mat = BlackBox_model.C
idmod.D_mat = BlackBox_model.D
idmod.L_mat = BlackBox_model.K
idmod.samp_T = 4 ;   % sampling time of 4 sec 

figure, plot( yk_id(:,1), 'k')
hold on
plot( yk(:,1),'xr')
hold off

figure, plot( yk_id(:,2), 'k')
hold on
hold on 
plot( yk(:,2),'xr')
hold off

uk = val_data.InputData ;
yk = val_data.OutputData ;

yk_val = idsim(BlackBox_model, uk);

figure, plot( yk_val(:,1), 'k')
hold on
plot( yk(:,1),'xr')
hold off

figure, plot( yk_val(:,2), 'k+')
hold on 
plot( yk(:,2),'xr')
hold off 

save_flag = input('Do you want to save the results? (Enter 1 for Yes and 0 for No) :') ;
if (save_flag ) 
     save TCL_BlackBox_OE_SS idmod 
end 


