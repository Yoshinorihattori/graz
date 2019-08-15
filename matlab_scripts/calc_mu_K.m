function mu = calc_mu_K( T )
    A = 1.1001E-4;
    B = 325.85;
    C = -207.30;

    mu = A .* exp( B ./ ( (T) + C) );
end
