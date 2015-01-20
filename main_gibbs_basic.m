function [phiChain, sig2Chain, tau2invChain] = main_gibbs_basic(data, lambda, nSim)

% Get the information of tha data 
[n,p] = size(data);                 % Retrieve the size of the data

% Allocate the matrix for simulation
phiChain   = zeros(p,p,nSim);       % phi-matrix x nSimulation
sig2Chain  = zeros(p,nSim);         % sig2       x nSimulation
tau2invChain = zeros(p,p,nSim);     % tau2inv    x nSimulation

% Deal with the first variable
yFirst      = data(:,1);
sig2First   = var(yFirst);
sig2Chain(1,:) = sig2Chain(1,:) + sig2First; 

% Initialize the starting values of the estimation

for k=2:p;
    
    [y,X] = getSubProblem(k,data);                  % k-th regression problem
 
    % Initialization of the starting points
    [ols_beta, ~, residual] = regress(y,X);
    
    phiChain(k,1:k-1,1) = -ols_beta';               % betas = OLS Beta 
    sig2Chain(k,1) = sumsqr(residual) / n;          % MLE variances for sigma
    tau2invChain(k,1:k-1,1) = 0.0001 * ones(1,k-1); % Use large number for tau2inv
    
    round = 2;
    
    while round ~= nSim
    
        % ================================= %
        %    START OF THE GIBBS SAMPLER
        % ================================= %
   
        % Do some preparation work here for easy coding
        % Retrieve the previous round of simulation results
        
        oldSig2 = sig2Chain(k,round-1);
        oldtau2Inv = tau2invChain(k,1:k-1,round-1);     % prevTau2Inv
        DINV = diag(oldtau2Inv);                        % diagnolized matrix
        
        try
            % Simulation for phi
                AINV = X'*X + DINV;
                RINV = chol(AINV);          
                R    = RINV\eye(k-1);
                m    = R * R' * X' * y;
                newPhi = m + R' * randn(k-1,1) * sqrt(oldSig2);

            % Simulation for sigma
                shape = (n+k-1)/2;
                scale = (sum((y-X*newPhi).^2) + newPhi' * DINV * newPhi)/2;
                newSig2 = 1 / gamrnd(shape, 1/scale);

            % Simulate tau2inverse
                lam_g = lambda^2 * ones(k-1,1);
                mu_g = sqrt(lambda^2 * newSig2 ./ newPhi.^2);
                newTau2inv = invGaussianMultipleParams(lam_g,mu_g,1);

            % save the new step
            phiChain(k,1:k-1,round) = newPhi';
            sig2Chain(k, round) = newSig2;
            tau2invChain(k,1:k-1,round) = newTau2inv;
            round = round + 1;
        catch
            round = round - 3;
        end
    end
    
end

end

