% Linear Quadratic Optimal Control Law implementation 

function [uk_LQR, LQR]  = LQOC_law( LQR, rk, yk, k ) 

     % Setpoint filtering 
     LQR.rk_f(:,k+1) = LQR.phy_r * LQR.rk_f(:,k) + (eye(LQR.n_op) - LQR.phy_r) * rk   ; 

     % Controller implementation at instant k 
     % Innovation/ MPM signal filtering 
     LQR.ek(:,k) = yk - LQR.C_mat * LQR.xkhat(:,k) ;  
     LQR.ek_f(:,k+1) = LQR.phy_e * LQR.ek_f(:,k) + (eye(LQR.n_op) - LQR.phy_e) * LQR.ek(:,k) ; 

     % Target State Computation 
     LQR.us_k(:,k) =  LQR.Gain_ss_inv * ( LQR.rk_f(:,k+1) - LQR.Gain_e * LQR.ek_f(:,k))  ;
     LQR.xs_k(:,k) = LQR.inv_mat * ( LQR.gam * LQR.us_k(:,k) + LQR.L_mat * LQR.ek_f(:,k))  ;

    % Control law computations 
    uk_LQR = LQR.us_k(:,k) - LQR.G_inf * (LQR.xkhat(:,k) - LQR.xs_k(:,k)) ;

    % Impose uk bound constraints
    uk_clipped = uk_bound_constraints( uk_LQR,  LQR.u_L, LQR.u_H ) ;
    uk_LQR = uk_clipped ;

    % Update the internal model state using the updated uk 
    LQR.xkhat(:,k+1) = LQR.phy * LQR.xkhat(:,k) + LQR.gam * uk_LQR + LQR.L_mat * LQR.ek(:,k) ;
end



