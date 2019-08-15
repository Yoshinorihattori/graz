function rho = calc_rho( T )
    A = 1045;
    B = -0.616;

    rho = A + B * (T+273.15);
end
