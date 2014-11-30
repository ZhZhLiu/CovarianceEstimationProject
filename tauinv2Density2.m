function [p] = tauinv2Density2(tau2invj, j, k, rho, phi, tau2inv, sig2,lam)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


if(j == 1)

   psi1 = phi(1) * sqrt(tau2invj);
   psi2 = phi(2) * sqrt(tau2inv(2));
   prod = psi1^2 - 2*rho * psi1 * psi2;
   
elseif(j == k-1)
   
    psikm1 = phi(k-1) * sqrt(tau2invj);
    psikm2 = phi(k-2) * sqrt(tau2inv(k-2));
    prod = psikm1^2 - 2*rho*psikm1 * psikm2;
    
else
    
    psiprev = phi(j-1) * sqrt(tau2inv(j-1));
    psij = phi(j) * sqrt(tau2invj);
    psinext = phi(j+1) * sqrt(tau2inv(j+1));
    prod = -2*rho*psiprev * psij + (1+rho^2) * psij^2 - 2*rho*psinext * psij;
    
end

    p = tau2invj^(-1.5) * exp( -prod / 2 / sig2 / (1-rho^2)) * exp( - lam^2 / 2 / tau2invj);
    


end

