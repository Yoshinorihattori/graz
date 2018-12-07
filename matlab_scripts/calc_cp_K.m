function cp = calc_cp_K( T )
    A = 0.818;
    B = 3.664E-3;

    cp = A + B * (T);
    cp =  cp *1000;
end
