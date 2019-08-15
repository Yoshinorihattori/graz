function cp = calc_cp( T )
    A = 0.818;
    B = 3.664E-3;

    cp = A + B * (T+273.15);
end
