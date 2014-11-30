function [sig2] = simulateSig2_rho(y, X, phi, CINV, n, k)


if(k==22)
    k;
end

shape = (n+k-1)/2;
scale = (sum((y-X*phi).^2) + phi' * CINV * phi)/2;
x = gamrnd(shape, 1/scale);
sig2 = 1/x;

end