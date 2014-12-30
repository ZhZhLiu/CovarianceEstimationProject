function [COV_EST, COVINV_EST, eLoss, qLoss] = main_rho_covEst(data, nSim, nBurnIn, lambda, SIG_TRUE, SIGINV)

    % RUN BASIC GIBBS SAMPLER
    [phiChain, sig2Chain, ~, rhoChain] = main_gibbs_rho(data, lambda, nSim);

    % PLOT SOME EVOLUTION GRAPH
    fileName = sprintf('Sig2Path_lambda_%4.2f', lambda);
    ExportSigma2Evolutiton(sig2Chain,fileName);
    
    fileName = sprintf('RhoPath_lambda_%4.2f', lambda);
    ExportRhoEvolutiton(rhoChain,fileName);

    % Compute the average of the results
    [COV_EST, COVINV_EST] = estimationFromGibbsPaths(phiChain, sig2Chain, nBurnIn);
    [eLoss,qLoss] = getLoss(COV_EST, SIG_TRUE);

    % Graph The Estimation
    fileName = sprintf('Result_lambda_%4.2f.eps',lambda);
    ExportCovarianceMatricesGraph(SIG_TRUE, COV_EST, SIGINV, COVINV_EST, data, fileName, lambda);

end

