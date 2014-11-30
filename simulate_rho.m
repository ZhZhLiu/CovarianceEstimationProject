function [ rho ] = simulate_rho(rho, phis, tau2invs, sig2s, p, dist)

 rhoCand = rho + dist * randn(1);
 
 if (rhoCand > 0 && rhoCand < 1)

     %desRatio = rho_density_ratio(rhoCand, rho, phis, tau2invs, sig2s, p);
     pdCand =  rho_density(rhoCand, phis, tau2invs, sig2s, p);
     pdOri =  rho_density(rho, phis, tau2invs, sig2s, p);

     acceptRatio = min ( 1,  pdCand / pdOri);
     
     if (rand(1) < acceptRatio)
         rho = rhoCand;
     end

 end

end

