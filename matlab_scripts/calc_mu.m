function mu = calc_mu( T )
    A = 1.1001E-4;
    B = 325.85;
    C = -207.30;

    mu = A .* exp( B ./ ( (T+273.15) + C) );
end
