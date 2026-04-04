% Function TCL_Dynamics.m
% Temperature Control Lab Experimental System Simulation 
% Function called by ODE solver for Dynamic Simulation 
% Given time and X(t), this function computes f(X(t), U(k), D(k)) 

function Zdot = TCL_JacobianFn( Z)
global TCL

TCL.T(1) = Z(1);  TCL.T(2 )= Z(2);
TCL.Uk(1) = Z(3) ; TCL.Uk(2) = Z(4) ; 

% Model differential equations 

QC12=TCL.U*TCL.As*(TCL.T(2)-TCL.T(1));
QR12=TCL.eps*TCL.sigma*TCL.A*(TCL.T(2)^4-TCL.T(1)^4);

Zdot(1) = ((TCL.U*TCL.A)/(TCL.m*TCL.Cp))*(TCL.Tinf-TCL.T(1)) ;
Zdot(1) = Zdot(1) + ((TCL.eps*TCL.sigma*TCL.A)/(TCL.m*TCL.Cp)) * (TCL.Tinf^4-TCL.T(1)^4) ;
Zdot(1) = Zdot(1) +TCL.alpha(1)*TCL.Uk(1) + QC12 + QR12;

Zdot(2) = ((TCL.U*TCL.A)/(TCL.m*TCL.Cp))*(TCL.Tinf-TCL.T(2)) ;
Zdot(2) = Zdot(2) + ((TCL.eps*TCL.sigma*TCL.A)/(TCL.m*TCL.Cp))*(TCL.Tinf^4-TCL.T(2)^4) ;
Zdot(2) = Zdot(2) + TCL.alpha(2)*TCL.Uk(2) - QC12 - QR12;

Zdot = Zdot';
end