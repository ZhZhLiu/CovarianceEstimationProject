function [newTau2inv] = simulateTau2Inv_rho(phi, sig2, rho, pretau2inv, k, lambda)

lam_g = lambda^2; 
newTau2inv = pretau2inv;   % Use the previous observation to serve as the default

for j = 1:k-1
        
    mu_g  = sqrt(lambda^2 * sig2 / phi(j)^2);       % Location parameter
    
    old  = pretau2inv(j);                           % Previous observation
    new  = ingrnd(mu_g,lam_g);                      % Candidate observation
    
    % This is the probability density
    p_old = tauinv2Density2(old, j, k, rho, phi, newTau2inv, sig2, lambda);
    p_new = tauinv2Density2(new, j, k, rho, phi, newTau2inv, sig2, lambda);
    
    % This is the transitional density
    K_old = invGauPdf(old,lam_g,mu_g);
    K_new = invGauPdf(new,lam_g,mu_g);
    
    % Get the acceptance probability
    acceptRatio = min ( 1,  p_old * K_new / p_new / K_old);
    
    % Update the chain the new position is accepted
    if (rand(1) < acceptRatio)
        newTau2inv(j) = new;
    end

end

end