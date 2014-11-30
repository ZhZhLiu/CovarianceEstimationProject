function [ ratio ] = rho_density_ratio(numRho, demRho, phis, tau2invs, sig2s, p)

psis = phis .* sqrt(tau2invs);

exponent = 0;

for k = 3:p 
    diff = psis(k,1:end-1) - numRho * psis(k,2:end);
    sss  = sum(diff.^2);
    exponentk = sss + (1-numRho^2) * psis(k,k-1)^2;
    exponent = exponent + exponentk / 2 / sig2s(k) / (1 - numRho^2);
end

exponent2 = 0;

for k = 3:p 
    diff = psis(k,1:end-1) - demRho * psis(k,2:end);
    sss  = sum(diff.^2);
    exponentk = sss + (1-demRho^2) * psis(k,k-1)^2;
    exponent2 = exponent2 + exponentk / 2 / sig2s(k) / (1 - demRho^2);
end

    
    ratio = ((1-numRho^2) / (1-demRho^2))^(- (p-2)*(p-1) / 4 ) * exp( - (exponent - exponent2));

end

