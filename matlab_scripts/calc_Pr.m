function Pr = calc_Pr( T )
    
    mu     = calc_mu( T );
    cp     = calc_cp( T );
    lambda = calc_lambda( T );

    Pr = ( mu .* cp .* 1000 ) ./ lambda;
end
