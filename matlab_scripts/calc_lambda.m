function lambda = calc_lambda( T )
    A = 0.157;
    B = -7.328E-5;

    lambda = A + B * (T+273.15);
end
