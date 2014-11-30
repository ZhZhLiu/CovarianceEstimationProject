%% Generate the dataset for estimation

clear

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

resultFileName = 'Loss_summary_2.txt';
myformat = '[Lam %5.3f]\t %7.4f\t %7.4f \n';
fid = fopen(resultFileName,'a');

n = 100;          % number of observations
p = 30;             % dimensions
sig2 = 0.01;        % variance of the data
rho = 0.9;          % exponential decay of the T matrix

% Use the producedure to generate the correct dataset c
[data, COVINV_TRUE, COV_TRUE, T_TRUE, DINV_TRUE] = generateDataset(n,p,sig2,rho,'decay');


%% Start some initialization process for the sampler

% Simulation Parameters
nSim = 1000;         % Number of Simulation
nBurnIn = 100;
lambda = 300;         % Regularization Parameter

% Allocate the matrix for simulation
phiChain   = zeros(p,p,nSim);       % phi-matrix x nSimulation
sig2Chain  = zeros(p,nSim);         % sig2       x nSimulation
rhoChain  = zeros(nSim,1);          % rho        x nSimulation
tau2invChain = zeros(p,p,nSim);     % tau2inv    x nSimulation

% Deal with the first variable
yFirst      = data(:,1);
sig2First   = var(yFirst);


% Initialize the starting values of the estimation
for k = 2:p
    
    % Use the regression results to initialize the Phi First
    [y,X] = getSubProblem(k,data);
    temp_beta =  regress(y,X)';
    phiChain(k,1:k-1,1) = temp_beta;
    
    % Use the the variance of the residual to initialize the sigma First
    sig2Chain(k,1) = sum(( y - X*temp_beta').^2) / n;
    
    % Use large number for tau2inv
    tau2invChain(k,1:k-1,1) = 0.0001 * ones(1,k-1);
    
    % Use the prior assumption of the rho to initialize the rho first
    rhoChain(1) = 0.2;
    
end


%% Starts the Gibbs Sampler for estimation

sigFirst = var(data(:,1));

for round = 1:nSim

    for k = 2:p
        
    % Get the hold of previous round of values
    preSig2 = sig2Chain(k,round);                   %Receive the previous round Sig2
    preTau2Inv = tau2invChain(k, 1:k-1, round);     %Receive the previous round Tau2Inv
    preRho = rhoChain(round);                     %Receive the previous round Rho

    % Get the sub problem for the k-th regression problem 
    [y,X] = getSubProblem(k,data);
        
    % Some transformation of the values, compute the CINV
    DINV_half = diag(preTau2Inv.^(0.5));    
    CINV = DINV_half * getLamMatrix(rho,k-1) * DINV_half;

    % Do the simulation
    newPhi  = simulatePhi_rho(y, X, sig2, CINV, k);         % Get the new observations for phi
    newSig2 = simulateSig2_rho(y, X, newPhi, CINV, n, k);   % Get the new observations for sig2
    
    if (k==2)
        newTau2Inv = simulateTau2Inv(newPhi, newSig2, k, lambda);
    else
        newTau2Inv = simulateTau2Inv_rho(newPhi, newSig2, preRho, preTau2Inv, k, lambda); % Simulate newTau2Inv
    end
    
    % Assign the new value to the chain
    phiChain(k,1:k-1,round+1) = newPhi';
    sig2Chain(k,round+1) = newSig2;
    tau2invChain(k,1:k-1,round+1) = newTau2Inv;
    
    end
    
    % Simulate rho
    dist = 0.1;
    rhoChain(round+1) = simulate_rho(preRho, phiChain(:,:,round+1),  tau2invChain(:,:,round+1), sig2Chain(:,round+1) , p, dist);
    %rhoChain(round+1) = rho;
end


%% Compute the average of the results

T_EST = -mean(phiChain(:,:,nBurnIn:end),3) + eye(p);
S_EST = median(sig2Chain(:,nBurnIn:end), 2);
S_EST(1) = sig2First;
DINV_EST = diag(1./S_EST);

COVINV_EST = T_EST' * DINV_EST * T_EST;
COV_EST = inv(COVINV_EST);
RHO_EST = mean(rhoChain(nBurnIn:end));

[eLoss,qLoss] = getLoss(COV_EST, COV_TRUE);

handle = figure(1);

subplot(2,3,1)
imagesc(COV_TRUE);
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
imagesc(COVINV_TRUE);
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


