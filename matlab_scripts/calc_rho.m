function rho = calc_rho( T )
    A = 1268.28;
    B = -0.66;

    rho = A + B * (T+273.15);
end
