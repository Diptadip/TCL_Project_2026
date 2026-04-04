% <---------------------------------------------------------------------------->
% function gradient.m
% This function computes numerical gradient matrix for any function vector  
% Input arguments to the function are
% fun_name  : String containing name of the MATLAB function,
%                     which returns function vector F(x) given x 
% x0             : vector at which gradient should be computed 
% Output arguments of the function are
% gradF        : (n x n) gradient matrix  
% ------------------------------------------------------------------------------

function gradF = Numerical_Jacobian( fun_name, Z0 ) 

nZ0 = length( Z0 ) ;

for i = 1 : nZ0 
    
    eps = Z0(i) / 100000;
    
    if ( eps == 0 )
        eps = 1 / 100000;
    end     
    
    Zp = Z0 ;
    Zp(i) = Zp(i) + eps ; 
    fp = feval( fun_name, Zp ) ;
    Zn = Z0 ;
    Zn(i) = Zn(i) - eps ; 
    fn = feval( fun_name, Zn ) ;   
    gradF(:,i) = (fp - fn) / (2 * eps) ;
    
end     
    