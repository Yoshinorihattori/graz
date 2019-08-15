function rho = calc_rho_K( T )
    A = 1268.28;
    B = -0.66;

    rho = A + B * (T);
end
