function [COV_EST, COVINV_EST, eLoss, qLoss, simChain] = main_basic_covEst(data, nSim, nBurnIn, lambda, SIG_TRUE, SIGINV, rootPath)

    % RUN BASIC GIBBS SAMPLER
    [phiChain, sig2Chain, tau2invChain] = main_gibbs_basic(data, lambda, nSim);

    % Compute the average of the results
    [COV_EST, COVINV_EST] = estimationFromGibbsPaths(phiChain, sig2Chain, nBurnIn);
    [eLoss,qLoss] = getLoss(COV_EST, SIG_TRUE);

    simChain.phi        = phiChain;
    simChain.sig2       = sig2Chain;
    simChain.tau2inv    = tau2invChain;
% 
%     % PLOT SOME EVOLUTION GRAPH
%     fileName = sprintf('Sig2Path_lambda_%4.2f', lambda);
%     ExportSigma2Evolutiton(sig2Chain,fileName, rootPath);
% 
%     % Graph The Estimation
%     fileName = sprintf('Result_lambda_%4.2f',lambda);
%     ExportCovarianceMatricesGraph(SIG_TRUE, COV_EST, SIGINV, COVINV_EST, data, lambda, fileName, rootPath);

end

