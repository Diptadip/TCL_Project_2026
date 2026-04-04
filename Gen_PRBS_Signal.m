function delta_uk = Gen_PRBS_Signal( n_ip, prbs_ampl )
% Function Gen_PRBS_Signal generated pseudo-random binary input sequence
%  Input arguments: current sample number., switching period, PRBS amplitude
% Output arguments: PRBS input

    ftmp = randn(n_ip,1) ;
    delta_uk = zeros(n_ip,1) ;
    for i = 1 : n_ip
        if ( ftmp(i) >0)
            delta_uk(i) = prbs_ampl(i) ;
        else
            delta_uk(i) = - prbs_ampl(i) ;
        end
    end
    

