
% Weighted least square objective function for estimation 
% parameters of nonlinear mechanistic model for TCL 
 
function obj_fn = obj_paraEstm_TCL(para)

 global TCL Yk_data Uk_data Xk_data
     
TCL.U = para(1);          % Overall Heat transfer coefficient
TCL.alpha(1) = para(2);  % Heater 1 efficiency parameter 
TCL.alpha(2) = para(3);  % Heater 2 efficiency parameter 

Xk(:,1) = Yk_data(:,1);  % Set initial state equal to the first set of measurements 

N_samp = length(Yk_data) ;   

obj_fn=0;

for k = 1:N_samp-1
    
    % Specify inputs at instant k 
    TCL.Uk(1)=Uk_data(1,k);  
    TCL.Uk(2)=Uk_data(2,k);

    % Process Simulation from sampling instant k to (k+1)
    
    % [t,Xt]=ode23('TCL_Dynamics',[0 TCL.Samp_T], Xk(:,k));  
    % Xk(:,k+1)=Xt(end,:)';
 
    % 2nd order RK method for numerical integration
    fk_0 = TCL_Dynamics( 0, Xk(:,k) ) ;
    Xk_1 = Xk(:,k) + TCL.Samp_T * fk_0 ; 
    fk_1 = TCL_Dynamics( 0, Xk_1 ) ;
    Xk(:,k+1) = Xk(:,k) + 0.5*TCL.Samp_T * (fk_0 + fk_1) ;

    Est_err_k = (Yk_data(:,k+1)-Xk(:,k+1))  ;
    obj_fn = obj_fn + Est_err_k'*inv(TCL.R)* Est_err_k ;

end

Xk_data = Xk ;
end