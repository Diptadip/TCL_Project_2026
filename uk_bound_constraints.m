% Function to impose input bounds

function uk_clipped = uk_bound_constraints( u_k, u_L, u_H )

% Ad-hoc implementation of bound constraints

uk_clipped = u_k ; 

if ( u_k(1) <= u_L(1) )   % Impose input constraints on uk(1)
    uk_clipped(1) = u_L(1) ;
elseif ( u_k(1) >= u_H(1) )
    uk_clipped(1) = u_H(1) ;
end
if ( u_k(2) <= u_L(2) ) % Impose output constraints on uk(2)
    uk_clipped(2) = u_L(2) ;
elseif ( u_k(2) >= u_H(2) )
    uk_clipped(2) = u_H(2) ;
end
end