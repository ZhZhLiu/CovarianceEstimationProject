function [newTau2inv] = simulateTau2Inv(phi, sig2, k, lambda)

newTau2inv = zeros(1,k-1);
for j = 1:k-1
    lam_g = lambda^2;
    mu_g  = sqrt(lambda^2 * sig2 / phi(j)^2);
    newTau2inv(1,j) = ingrnd(mu_g,lam_g);
end

end