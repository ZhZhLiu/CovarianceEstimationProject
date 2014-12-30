%% Generate the dataset for estimation

clear

s = RandStream('mt19937ar','Seed',1035);
RandStream.setGlobalStream(s);

resultFileName = 'Loss_summary_2.txt';
myformat = '[Lam %5.3f]\t %7.4f\t %7.4f \n';
fid = fopen(resultFileName,'a');

n = 100;          % number of observations
p = 30;             % dimensions
sig2 = 0.01;        % variance of the data
rho = 0.5;          % exponential decay of the T matrix

% Use the producedure to generate the correct dataset c
[data, COVINV_TRUE, COV_TURE, T_TRUE, DINV_TRUE] = generateDataset(n,p,sig2,rho,'decay');

%% Start the Gibbs Sampler for estimation procedure

% Simulation Parameters
nSim = 5000;         % Number of Simulation
nBurnIn = 1000;      % Burn in number

lambdas = linspace(0,20,101);
lambdas = lambdas(2:end);

allLosses = zeros(2,length(lambdas));

for lambda = lambdas
    i = 1;
    disp(sprintf('lambda = %4.2f',lambda));
    [COV_EST, COVINV_EST, eLoss,qLoss ] = gibbsSample_Fundamental(data, nSim, nBurnIn, lambda, COV_TURE, COVINV_TRUE);
    fprintf(fid, myformat, lambda, eLoss, qLoss);
    allLosses(:,i) = [eLoss;qLoss];
    i=i+1;
end

fclose(fid);