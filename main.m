%% Generate the dataset for estimation

clear

n = 100;          % number of observations
p = 30;             % dimensions
sig2 = 0.01;        % variance of the data
rho = 0.8;          % exponential decay of the T matrix

% Use the producedure to generate the correct dataset c
[data,SIGINV, SIG_TRUE, T_TRUE, DINV_TRUE] = generateDataset(n,p,sig2,rho,'cutoff');

%% Start the Gibbs Sampler for estimation procedure

% Simulation Parameters
nSim = 5000;         % Number of Simulation
nBurnIn = 1000;
lambda = 2;         % Regularization Parameter

% Allocate the matrix for simulation
phiChain   = zeros(p,p,nSim);       % phi-matrix x nSimulation
sig2Chain  = zeros(p,nSim);         % sig2       x nSimulation
tau2invChain = zeros(p,p,nSim);     % tau2inv    x nSimulation

% Deal with the first variable
yFirst      = data(:,1);
sig2First   = var(yFirst);

% Initialize the starting values of the estimation
for k = 2:p
    
    % Use the regression results to initialize the Phi First
    [y,X] = getSubProblem(k,data);
    [ols_beta, ~, residual] = regress(y,X);
    phiChain(k,1:k-1,1) = -ols_beta';
    
    % Use the the variance of the residual to initialize the sigma First
    sig2Chain(k,1) = sumsqr(residual) / n;
 
    % Use large number for tau2inv
    tau2invChain(k,1:k-1,1) = 0.0001 * ones(1,k-1);
end


for k=2:p;
    
    [y,X] = getSubProblem(k,data);                  % k-th regression problem
 
    for round = 2:nSim;
    
        % Do some preparation work here for easy coding
        % Retrieve the previous round of simulation results
        
        oldSig2 = sig2Chain(k,round-1);                 % prevSig2 
        oldtau2Inv = tau2invChain(k,1:k-1,round-1);     % prevTau2Inv
        DINV = diag(oldtau2Inv);                        % diagnolized matrix
        
        % Start the gibbs sampler
        
        % Simulation for phi
            AINV = X'*X + DINV;
            RINV = chol(AINV);
            R    = RINV\eye(k-1);
            m    = R * R' * X' * y;
            newPhi = m + R' * randn(k-1,1) * sqrt(sig2);

        % Simulation for sigma
            shape = (n+k-1)/2;
            scale = (sum((y-X*newPhi).^2) + newPhi' * DINV * newPhi)/2;
            newSig2 = 1 / gamrnd(shape, 1/scale);
        
        % Simulate tauinver
            
            newTau2inv = zeros(1,k-1);
            for j = 1:k-1
                lam_g = lambda^2;
                mu_g  = sqrt(lambda^2 * sig2 / phi(j)^2);
                newTau2inv(1,j) = ingrnd(mu_g,lam_g);
            end
        
        
        newTau2inv = simulateTau2Inv(newPhi, newSig2, k, lambda);   % Simulation of tau2Inv

        % save the new step
        phiChain(k,1:k-1,round) = newPhi';
        sig2Chain(k, round) = newSig2;
        tau2invChain(k,1:k-1,round) = newTau2inv;

    end
    
end

%% Compute the average of the results

T_EST = -mean(phiChain(:,:,nBurnIn:end),3) + eye(p);
S_EST = mean(sig2Chain(:,nBurnIn:end), 2);
S_EST(1) = sig2First;
DINV_EST = diag(1./S_EST);

COVINV_EST = T_EST' * DINV_EST * T_EST;
COV_EST = inv(COVINV_EST);

%% COMPUTE THE LOSS

[eLoss,qLoss] = getLoss(COV_EST, SIG_TRUE);

handle = figure(1);

subplot(2,3,1)
imagesc(SIG_TRUE);
title('True Covariance')
axis square

subplot(2,3,2)
imagesc(cov(data));
title('Sample Covariance')
axis square

subplot(2,3,3)
imagesc(COV_EST);
title(sprintf('Regularized Covariance \\lambda=%4.2f', lambda))
axis square

subplot(2,3,4)
imagesc(SIGINV);
title('True Covariance Inverse')
axis square

subplot(2,3,5)
imagesc(inv(cov(data)));
title('Sample Covariance Inverse')
axis square

subplot(2,3,6)
imagesc(COVINV_EST);
title(sprintf('Regularized Covariance Inverse \\lambda=%4.2f', lambda))
axis square

nameFile_EPS = sprintf('Figures/Result_lambda_%4.2f.eps', lambda);
nameFile_FIG = sprintf('Figures/Result_lambda_%4.2f.fig', lambda);


set(gcf,'PaperPositionMode','auto')
set(handle, 'Position', [100, 100, 1400, 900]);
%print('-dpsc2','-zbuffer','-r200')
print(handle, '-depsc', nameFile_EPS);
savefig(handle,nameFile_FIG);

