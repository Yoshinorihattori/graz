function lambda = calc_lambda( T )
    A = 0.2134;
    B = 6.071E-4;

    lambda = A + B * (T+273.15);
end
