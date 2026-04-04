%TC LAB open loop simulation and comparison with experiment data

clc, clear all, close all

global TCL

load TCL_MechModel_Parameters     % Load TCL Nonlinear Model parameters 

n_ip = 2 ; n_op = 2 ; , n_st = 2 ;

% Find Equilibrium/ Steady State Operating Point for the specified inputs 

Xs0 = TCL.Xs ;       % Initial Guess
oldopts = optimset ;
tolerance = 1e-10 ;
options = optimset(oldopts,'MaxIter', 1e6, 'Display','off', 'TolFun',tolerance, 'TolX', tolerance) ;
TCL.Xs = fsolve('TCL_SteadyState', Xs0, options ) ;

fprintf('\n Steady State corresponding to the specified input is as follows: ') ;
display(TCL.Xs')
% fprintf('\n Hit ANY Key to Continue... \n') ;pause

% Generate linear perturbation model at the current steady state 
Z_vec = [ TCL.Xs' TCL.Us' ]' ;
Jacob_mat = Numerical_Jacobian( 'TCL_JacobianFn', Z_vec ) ;
% Continuous time linear perturbation model 
TCL.A_mat = Jacob_mat(:,1:n_st) ; 
TCL.B_mat  = Jacob_mat(:,n_st+1:n_st+n_ip) ;

% Create a state space object in Matlab 
cmod = ss( TCL.A_mat, TCL.B_mat, TCL.C_mat, TCL.D_mat ) ;     

fprintf('\n Continuous Time Linear State Space Model at chosen equilibrium point') ;
cmod
% fprintf('\n Hit ANY Key to Continue... \n') ;pause

%  Generate discrete time linear perturbation model 
dmod = c2d( cmod, TCL.Samp_T ) ;          % Continuous to discrete conversion 
fprintf('\n Discrete Time Linear State Space Model at the chosen equilibrium point') ;
dmod
% fprintf('\n Hit ANY Key to Continue... \n') ;pause

TCL.phy = dmod.A ; 
TCL.gama = dmod.B ;   

save TCL_MechModel_Parameters TCL

