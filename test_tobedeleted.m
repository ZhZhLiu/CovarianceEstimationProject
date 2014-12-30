%% RANDOMIZATION OF SOME MATRIX
clear;

oldRho    = rand();
canRho = rand();
sig2 = 0.1;

k = 10;



LAMBDA_UPPER = - diag(ones(k-2,1) * oldRho, 1);
LAMBDA_LOWER = - diag(ones(k-2,1) * oldRho, -1);
LAMBDA_OLD = eye(k-1) * (1 + oldRho^2) + LAMBDA_UPPER + LAMBDA_LOWER;
LAMBDA_OLD(1,1) = 1;
LAMBDA_OLD(end,end) = 1;
LAMBDA_OLD = LAMBDA_OLD / ( 1 - oldRho^2);

LAMBDA_UPPER = - diag(ones(k-2,1) * canRho, 1);
LAMBDA_LOWER = - diag(ones(k-2,1) * canRho, -1);
LAMBDA_CAN = eye(k-1) * (1 + canRho^2) + LAMBDA_UPPER + LAMBDA_LOWER;
LAMBDA_CAN(1,1) = 1;
LAMBDA_CAN(end,end) = 1;
LAMBDA_CAN = LAMBDA_CAN / ( 1 - canRho^2);


phi         = rand(k-1, 1) * 0.2 + 1;
tau2Inv     = rand(k-1, 1) * 0.2 + 1;
DINV_Half   = sqrt(diag(tau2Inv));

CINV_OLD = DINV_Half * LAMBDA_OLD * DINV_Half;
CINV_CAN = DINV_Half * LAMBDA_CAN * DINV_Half;

DET_C_OLD = (1-oldRho^2)^(k-2) / prod(tau2Inv);
DET_C_CAN = (1-canRho^2)^(k-2) / prod(tau2Inv);

density_oldRho = 1 / sqrt(DET_C_OLD) * exp(- phi' * CINV_OLD * phi / 2 / sig2 );
density_canRho = 1 / sqrt(DET_C_CAN) * exp(- phi' * CINV_CAN * phi / 2 / sig2 );
