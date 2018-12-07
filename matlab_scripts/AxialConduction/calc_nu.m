function nu = calc_nu( T )
    rho=calc_rho(T);
    mu=calc_mu(T);
    nu = mu./rho;
end
