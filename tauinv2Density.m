function [pp] = tauinv2Density(tau2cand,j, phi, tau2inv, sig2, rho, lambda, k)

if( j == 1)
    
    tauinv2 = tau2cand;

    psi1 = phi(1) * sqrt(tau2cand);
    psi2 = phi(2) * sqrt(tau2inv(2));
    
    exponent = exp ( - (psi1 - rho * psi2)^2 / 2 / sig2 / (1-rho^2));
    expdes   = exp ( - lambda^2 / 2 / tauinv2);
    
    praw =  (tauinv2)^0.5 * exponent * expdes;
   
elseif(j == k-1)

    tauinv2 = tau2cand;

    psik1 = phi(j) * sqrt(tau2cand);
    psik2 = phi(j-1) * sqrt(tau2inv(j-1));

    exponent = exp ( - (psik1 - rho * psik2)^2 / 2 / sig2 / (1-rho^2));
    expdes   = exp ( - lambda^2 / 2 / tauinv2);
    
    praw =  (tauinv2)^0.5 * exponent * expdes;

else
    
    tauinv2 = tau2cand;
    
    psi_prev = phi(j-1) * sqrt(tau2inv(j-1));
    psi = phi(j) * sqrt(tau2cand);
    psi_next = phi(j+1) * sqrt(tau2inv(j+1));
    
    lagDiff1 = psi_prev - rho * psi;
    lagDiff2 = psi - rho* psi_next;

    
    exponent = exp ( - (lagDiff1^2 + lagDiff2^2) / 2 / sig2 / (1-rho^2) );
    expdes   = exp ( - lambda^2 / 2 / tauinv2);
    praw =  (tauinv2)^0.5 * exponent * expdes;
    
end

    pp = praw * tau2cand^(-2);

end

