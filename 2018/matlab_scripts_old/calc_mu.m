function mu = calc_mu( T )
    A = 5.894E-5;
    B = 857.4;
    C = -172.2;

    mu = A .* exp( B ./ ( (T+273.15) + C) );
end
