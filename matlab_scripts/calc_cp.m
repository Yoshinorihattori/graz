function cp = calc_cp( T )
    A = 2.0148;
    B = 4.50E-3;

    cp = A + B * (T+273.15);
end
