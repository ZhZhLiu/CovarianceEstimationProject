function [sig2] = simulateSig2(y, X, phi, tau2Inv, n, k)

DINV = diag(tau2Inv);
shape = (n+k-1)/2;
scale = (sum((y-X*phi).^2) + phi' * DINV * phi)/2;
x = gamrnd(shape, 1/scale);
sig2 = 1/x;

end