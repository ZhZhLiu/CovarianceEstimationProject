function [COV_EST, COVINV_EST, eLoss,qLoss ] = gibbsSample_Fundamental(data, nSim, nBurnIn, lambda, COV_TRUE, COVINV_TRUE)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[n,p] = size(data);

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
    temp_beta =  regress(y,X)';
    phiChain(k,1:k-1,1) = -temp_beta;
    
    % Use the the variance of the residual to initialize the sigma First
    sig2Chain(k,1) = sum(( y - X*temp_beta').^2) / n;
    
    % Use large number for tau2inv
    tau2invChain(k,1:k-1,1) = 0.0001 * ones(1,k-1);

end


for k=2:p;
    display(k);
    [y,X] = getSubProblem(k,data);                  % k-th regression problem
 
    for round = 2:nSim;
    
        % Do some preparation work here for easy coding

        oldSig2 = sig2Chain(k,round-1);                 % prevSig2 
        oldtau2Inv = tau2invChain(k,1:k-1,round-1);     % prevTau2Inv

        % Start simulation one by one
        newPhi = simulatePhi(y,X,oldSig2, oldtau2Inv, k);        % Simulation of phi            
        newSig2= simulateSig2(y,X, newPhi, oldtau2Inv, n, k);       % Simulation of sig2
        newTau2inv = simulateTau2Inv(newPhi, newSig2, k, lambda);   % Simulation of tau2Inv

        % save the new step
        phiChain(k,1:k-1,round) = newPhi';
        sig2Chain(k, round) = newSig2;
        tau2invChain(k,1:k-1,round) = newTau2inv;

        
        
    end
    
    
        handle = figure(1);
        if(k>0)    
            
            subplot(3,1,1)
            plot(sig2Chain(k,1:round));
            xlim([1,nSim]);
            title('MC Paths of \sigma^2');
            
            subplot(3,1,2)
            rPhi =  reshape(phiChain(k,1:k-1,:),[k-1,nSim]);
            plot(rPhi(:,1:round)');
            xlim([1,nSim]);
            title('MC Paths of \phi');
            
            subplot(3,1,3)
            rPhi =  reshape(tau2invChain(k,1:k-1,:),[k-1,nSim]);
            plot(rPhi(:,1:round)');
            xlim([1,nSim]);
            title('MC Paths of 1/\tau^2');
            
        end
    
    
    nameFile_EPS = sprintf('testConvergence_rho%4.2f_k%d.eps',0.7,k);
    nameFile_FIG = sprintf('testConvergence_rho%4.2f_k%d.jpg',0.7,k);
    print(handle, '-depsc', nameFile_EPS);
    print(handle, '-djpeg', nameFile_FIG);
    
end


%% Compute the average of the results

T_EST = -mean(phiChain(:,:,nBurnIn:end),3) + eye(p);
S_EST = mean(sig2Chain(:,nBurnIn:end), 2);
S_EST(1) = sig2First;
DINV_EST = diag(1./S_EST);

COVINV_EST = T_EST' * DINV_EST * T_EST;
COV_EST = inv(COVINV_EST);

%% COMPUTE THE LOSS

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

end

