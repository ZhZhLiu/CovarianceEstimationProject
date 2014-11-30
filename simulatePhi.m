function [newPhi] = simulatePhi(y, X, sig2, tau2Inv, k)

DINV = diag(tau2Inv);

AINV = X'*X + DINV;
RINV = chol(AINV);
R    = RINV\eye(k-1);
m    = R * R' * X' * y;

newPhi = m + R' * randn(k-1,1) * sqrt(sig2);

end