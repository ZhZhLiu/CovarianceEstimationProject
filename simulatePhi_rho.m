function [newPhi] = simulatePhi_rho(y, X, sig2, CINV, k)

AINV = X'*X + CINV;
RINV = chol(AINV);
R    = RINV\eye(k-1);
m    = R * R' * X' * y;
newPhi = m + R' * randn(k-1,1) * sqrt(sig2);

end