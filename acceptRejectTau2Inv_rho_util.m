function [ updatedTau2Inv ] = acceptRejectTau2Inv_rho_util(tau2Invs, mus, lams, rho, sig2, phi)

[~, regNum] = size(tau2Invs);

acceptRand = rand(1,regNum);        % Uniformly distributed random variable for acceptance rejection
phi = repmat(phi,2,1);

kernel_density = invGaussianPdfMultipleParams(tau2Invs, mus, lams);
taus = sqrt(1./tau2Invs);
psi  = sqrt(phi ./ taus);

for j = 1:regNum

    tau_j = taus(:,j);
    kernel = kernel_density(:,j);
    
    if j == 1
        inner = (psi(:,j) - rho * psi(1,j+1)).^2;
    elseif j == regNum
        inner = (1-rho^2) * psi(:,j).^2;
    else
        inner = (psi(:,j) - rho * psi(1,j+1)).^2 + (psi(1,j-1) - rho * psi(:,j)).^2;
    end

    exponent = - inner/ ( 2 * sig2* (1-rho^2)) - lams(j)^2 / 2 * tau_j.^2;
    
    p = exp(exponent) .* tau_j.^(-2.5);
    alpha = min(1, p(1) * kernel(2) / p(2) / kernel(1));

    if(acceptRand(j) < alpha) 
        tau2Invs(1,j) = tau2Invs(2,j);
        psi(1,j) = psi(2,j);
    end
end

updatedTau2Inv = tau2Invs(1,:);

end

