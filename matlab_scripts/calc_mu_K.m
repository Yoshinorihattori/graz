function mu = calc_mu_K( T )
    A = 5.894E-5;
    B = 857.4;
    C = -172.2;

    mu = A .* exp( B ./ ( (T) + C) );
end
