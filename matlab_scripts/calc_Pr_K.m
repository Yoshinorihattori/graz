function Pr = calc_Pr_K( T )
    
    mu     = calc_mu( T );
    cp     = calc_cp( T );
    lambda = calc_lambda( T );

    Pr = ( mu .* cp ) ./ lambda;
end
