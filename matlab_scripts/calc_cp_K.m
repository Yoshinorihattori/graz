function cp = calc_cp_K( T )
    A = 2.0148;
    B = 4.50E-3;

    cp = A + B * (T);
    cp =  cp *1000;
end
