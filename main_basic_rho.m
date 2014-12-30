%% Generate the dataset for estimation

clear

n = 100;            % number of observations
p = 30;             % dimensions
sig2 = 0.01;        % variance of the data
rho = 0.5;          % exponential decay of the T matrix

% Use the producedure to generate the correct dataset c
[data,SIGINV, SIG_TRUE, T_TRUE, DINV_TRUE] = generateDataset(n,p,sig2,rho,'decay');

%% Start the Gibbs Sampler for estimation procedure

% Simulation Parameters
nSim = 5000;        % Number of Simulation
nBurnIn = 500;
lambda = 20;         % Regularization Parameter

%% RUN BASIC GIBBS SAMPLER

[phiChain, sig2Chain, tau2invChain, rhoChain] = main_gibbs_rho(data, lambda, nSim);

%% PLOT SOME EVOLUTION GRAPH
fileName = sprintf('Sig2Path_lambda_%4.2f', lambda);
ExportSigma2Evolutiton(sig2Chain,fileName);

fileName = sprintf('RhoPath_lambda_%4.2f', lambda);
ExportRhoEvolutiton(rhoChain,fileName);

%% Compute the average of the results

[COV_EST, COVINV_EST] = estimationFromGibbsPaths(phiChain, sig2Chain, nBurnIn);
[eLoss,qLoss] = getLoss(COV_EST, SIG_TRUE);


%% Graph The Estimation

fileName = sprintf('Result_lambda_%4.2f.eps',lambda);
ExportCovarianceMatricesGraph(SIG_TRUE, COV_EST, SIGINV, COVINV_EST, data, fileName, lambda);
