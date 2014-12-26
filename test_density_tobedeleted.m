clear;

k = 11;

rho = 0.6;
sig2 = 0.1;
lambda = 2;

acceptRand = rand(1,k-1);

phi = 0.8 + rand(1,k-1) * 0.15;

mus         = 1 + rand(1,k-1) * 0.1;
lams        = 1 + rand(1,k-1) * 0.1;

oldTau2Inv  = 1 + rand(1,k-1) * 0.15;
canTau2Inv  = 1 + rand(1,k-1) * 0.15;

tau2Inv = [oldTau2Inv; canTau2Inv];


acceptRejectTau2Inv_rho_util(tau2Inv, mus, lams, rho, sig2, phi)


%%

j = 5;

for j = 1:k-1

    tau_j = tau(:,j);
    
    if j == 1
        inner = (psi(:,j) - oldRho * psi(1,j+1)).^2;
    elseif j == k-1
        inner = (1-oldRho^2) * psi(:,j).^2;
    else
        inner = (psi(:,j) - oldRho * psi(1,j+1)).^2 + (psi(1,j-1) - oldRho * psi(:,j)).^2;
    end

    exponent = - inner/ ( 2 * sig2* (1-oldRho^2)) - lambda^2 / 2 * tau_j.^2;
    density = exp(exponent) .* tau_j.^(-2.5);
    kernel = K(:,j);

    alpha = min(1, density(1) * kernel(2) / density(2) / kernel(1));
    a = acceptRand(j);
    
    if(a < alpha) 
        tau2Inv(1,j) = tau2Inv(2,j);
        psi(1,j) = psi(2,j);
    end
end

tau2Inv