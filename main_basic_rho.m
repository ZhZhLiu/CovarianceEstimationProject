clear

%% Generate the dataset for estimation

% The following are the parameters for running the simulation

n       = 100;          % number of observations
p       = 30;           % dimensions
sig2    = 0.01;         % variance of the data
rho     = 0.5;          % exponential decay of the T matrix

% Use the parameter to generate the random dataset 
datasetType = 'decay';
[data,SIGINV, SIG_TRUE, T_TRUE, DINV_TRUE] = generateDataset(n,p,sig2,rho,'datasetType');

% Parameters for gibbs sampler

nSim = 500;        % Number of simulation
nBurnIn = 50;      % Number of the simulated observations to be discarded
lambda = 1;        % Regularization Parameter

% Define the path for exporting the results
SIMULATION_RESULT_PATH = sprintf('Results/Rho/%s_%4.2f_%4.2f_%d_%d/', datasetType, rho, sig2, n, p);
if exist(SIMULATION_RESULT_PATH,'dir') == 0
    mkdir(SIMULATION_RESULT_PATH);
end

%% RUN BASIC GIBBS SAMPLER
[phiChain, sig2Chain, tau2invChain, rhoChain] = main_gibbs_rho(data, lambda, nSim);

%% PLOT SOME EVOLUTION GRAPHS
fileNamePattern = sprintf('Sig2Path_lambda_%4.2f', lambda);
ExportSigma2Evolutiton(sig2Chain, fileNamePattern, SIMULATION_RESULT_PATH);

fileName = sprintf('RhoPath_lambda_%4.2f', lambda);
ExportRhoEvolutiton(rhoChain,fileNamePattern, SIMULATION_RESULT_PATH);

%% Compute the average of the results

[COV_EST, COVINV_EST] = estimationFromGibbsPaths(phiChain, sig2Chain, nBurnIn);
[eLoss,qLoss] = getLoss(COV_EST, SIG_TRUE);

%% Graph The Estimation
fileNamePattern = sprintf('Result_lambda_%4.2f',lambda);
ExportCovarianceMatricesGraph(SIG_TRUE, COV_EST, SIGINV, COVINV_EST, data, lambda, fileNamePattern, SIMULATION_RESULT_PATH);
