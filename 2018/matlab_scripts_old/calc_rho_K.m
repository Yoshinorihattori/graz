function rho = calc_rho_K( T )
    A = 1045;
    B = -0.616;

    rho = A + B * (T);
end
