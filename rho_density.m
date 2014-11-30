function [ des ] = rho_density(rho, phis, tau2invs, sig2s, p)

psis = phis .* sqrt(tau2invs);

exponent = 0;

for k = 3:p 
    diff = psis(k,1:end-1) - rho * psis(k,2:end);
    sss  = sum(diff.^2);
    exponentk = sss + (1-rho^2) * psis(k,k-1)^2;
    exponent = exponent + exponentk / 2 / sig2s(k) / (1 - rho^2);
end
    
    des = (1-rho^2)^(- (p-2)*(p-1) / 4 ) * exp( - exponent);

end

