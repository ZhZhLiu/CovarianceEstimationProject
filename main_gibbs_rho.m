function [phiChain, sig2Chain, tau2invChain, rhoChain] = main_gibbs_rho(data, lambda, nSim)

% Get the information of tha data 
[n,p] = size(data);                 % Retrieve the size of the data

% Allocate the matrix for simulation
rhoChain   = zeros(1,nSim);         % rho        x nSimulation
phiChain   = zeros(p,p,nSim);       % phi-matrix x nSimulation
sig2Chain  = zeros(p,nSim);         % sig2       x nSimulation
tau2invChain = zeros(p,p,nSim);     % tau2inv    x nSimulation

% Deal with the first variable
yFirst      = data(:,1);
sig2First   = var(yFirst);
sig2Chain(1,:) = sig2Chain(1,:) + sig2First; 

% Initialize the starting values of the estimation
for k=2:p
    
    [y,X] = getSubProblem(k,data);                  % k-th regression problem
 
    % Initialization of the starting points
    [ols_beta, ~, residual] = regress(y,X);
    
    phiChain(k,1:k-1,1) = -ols_beta';               % betas = OLS Beta 
    sig2Chain(k,1) = sumsqr(residual) / n;          % MLE variances for sigma
    tau2invChain(k,1:k-1,1) = 0.0001 * ones(1,k-1); % Use large number for tau2inv
end


rhoRand = randn(1,nSim);

for round = 2:nSim;

    oldRho = rhoChain(round-1);
    canRho = oldRho + 0.05 * rhoRand(round);
    canRho = canRho * (abs(canRho) < 1) + oldRho * (abs(canRho) >= 1);
    
    logDensityRatio = 0;

    for k=2:p;
    
        [y,X] = getSubProblem(k,data);
        
        % ================================= %
        %    START OF THE GIBBS SAMPLER
        % ================================= %
   
        % Initialize the lambda inverse matrix       
        LINV_Diag = ones(k-1,1) * (1+oldRho^2);
        LINV_Diag(1) = 1;
        LINV_Diag(end) = 1;
        LINV = diag(ones(k-2,1) * -oldRho, 1) + diag(ones(k-2,1) * -oldRho, -1) +  diag(LINV_Diag);
        LINV = LINV / (1-oldRho^2);
    
        
        % Do some preparation work here for easy coding
        % Retrieve the previous round of simulation results
        oldSig2     = sig2Chain(k,round-1);
        oldTau2Inv  = tau2invChain(k,1:k-1,round-1);      % prevTau2Inv
        DINV        = diag(oldTau2Inv);                   % diagnolized matrix
        DINV_Half   = sqrt(DINV);                         % square-root of DINV
        CINV        = DINV_Half * LINV * DINV_Half;       % imposed correlation
        
        
        % Simulation for phi
            AINV = X'*X + CINV;
            RINV = chol(AINV);
            R    = RINV\eye(k-1);
            m    = R * R' * X' * y;
            newPhi = m + R' * randn(k-1,1) * sqrt(oldSig2);

        % Simulation for sigma
            shape = (n+k-1)/2;
            scale = (sum((y-X*newPhi).^2) + newPhi' * CINV * newPhi)/2;
            newSig2 = 1 / gamrnd(shape, 1/scale);
        
        % Simulate tau2inverse
            lam_g = lambda^2 * ones(k-1,1);
            mu_g = sqrt(lambda^2 * newSig2 ./ newPhi.^2);
            canTau2Inv = invGaussianMultipleParams(lam_g,mu_g,1);
            tau2Inv = [oldTau2Inv; canTau2Inv];
            
            if k==2
                newTau2Inv = canTau2Inv;
            else
                newTau2Inv = acceptRejectTau2Inv_rho_util(tau2Inv, mu_g, lam_g, oldRho, newSig2, newPhi');
            end
            
        % save the new step
        phiChain(k,1:k-1,round) = newPhi';
        sig2Chain(k, round) = newSig2;
        tau2invChain(k,1:k-1,round) = newTau2Inv;
        
        % calculate the density of rho
        
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

        DINV_Half = sqrt(diag(newTau2Inv));

        CINV_OLD = DINV_Half * LAMBDA_OLD * DINV_Half;
        CINV_CAN = DINV_Half * LAMBDA_CAN * DINV_Half;

        DET_C_OLD = (1-oldRho^2)^(k-2) / prod(newTau2Inv);
        DET_C_CAN = (1-canRho^2)^(k-2) / prod(newTau2Inv);

 %       density_oldRho = density_oldRho  / sqrt(DET_C_OLD) * exp(- newPhi' * CINV_OLD * newPhi / 2 / newSig2 );
 %       density_canRho = density_canRho  / sqrt(DET_C_CAN) * exp(- newPhi' * CINV_CAN * newPhi / 2 / newSig2 );
    
        logDensityRatio_this = -newPhi'*(CINV_CAN - CINV_OLD) * newPhi / 2 / newSig2 - 0.5 * log(DET_C_CAN / DET_C_OLD);
        logDensityRatio = logDensityRatio + logDensityRatio_this;
    end
   
    alphaTest = min(1, exp(logDensityRatio));
   
    if rand < alphaTest
        rhoChain(round) = canRho;
    else
        rhoChain(round) = oldRho;
    end
%     s = sprintf('OLDRHO:%6.4f  CANRHO:%6.4f, Accept:%6.4f, alpha:%10.5f ', oldRho, canRho, rhoChain(round), alphaTest);
%     disp(s)
end

end

