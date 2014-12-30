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
nSim = 5000;        % Number of Simulation;
nBurnIn = 100;      % Number of Observations to burn-in;

lambdas = linspace(0,10,101);
lambdas = lambdas(2:end);
nLam = length(lambdas);

%% RUN BASIC GIBBS SAMPLER

% Initialize the Results

COV_EST_RESULT = zeros(p,p,nLam);
COVINV_EST_RESULT = zeros(p,p,nLam);
eLosses = zeros(nLam,1);
qLosses = zeros(nLam,1);

parfor i = 1:nLam
    [COV_EST(:,:,i) , COVINV_EST(:,:,i), eLosses(i), qLoss(i)] = main_rho_covEst(data, nSim, nBurnIn, lambdas(i), SIG_TRUE, SIGINV);
end

%% SAVE WORKSPACE
save('SaveWorkSpace/dataRun.mat');

