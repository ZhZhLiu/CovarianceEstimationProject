%% Generate the dataset for estimation
clear

% The following are the parameters for running the simulation

n       = 100;          % number of observations
p       = 30;           % dimensions
sig2    = 0.01;         % variance of the data
rho     = 0.5;          % exponential decay of the T matrix

% Use the parameter to generate the random dataset 
datasetType = 'decay';
[data,SIGINV, SIG_TRUE, T_TRUE, DINV_TRUE] = generateDataset(n,p,sig2,rho,'datasetType');

% Define the path for exporting the results
SIMULATION_RESULT_PATH = sprintf('Results/Basic/%s_%4.2f_%4.2f_%d_%d/', datasetType, rho, sig2, n, p);
if exist(SIMULATION_RESULT_PATH,'dir') == 0
    mkdir(SIMULATION_RESULT_PATH);
end

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
    [COV_EST(:,:,i) , COVINV_EST(:,:,i), eLosses(i), qLoss(i)] = main_rho_covEst(data, nSim, nBurnIn, lambdas(i), SIG_TRUE, SIGINV, SIMULATION_RESULT_PATH);
end

%% SAVE WORKSPACE
WORKSPACE_LOCATION_PATH = strcat(SIMULATION_RESULT_PATH, 'WORKSPACE/');
if exist(WORKSPACE_LOCATION_PATH, 'dir') == 0
    mkdir(WORKSPACE_LOCATION_PATH)
end
save(strcat(WORKSPACE_LOCATION_PATH,'dataRun.mat'));

%% PLOT LOSS FILES
handle = figure(5);
subplot(1,2,1);
title('Entropy Losses');
plot(lambdas,eLosses);
xlabel('\lambda');
subplot(1,2,2);
title('Quadratic Losses');
plot(lambdas,qLosses);
xlabel('\lambda');