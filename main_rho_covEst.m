function [COV_EST, COVINV_EST, eLoss, qLoss] = main_rho_covEst(data, nSim, nBurnIn, lambda, SIG_TRUE, SIGINV, rootPath)

    % RUN BASIC GIBBS SAMPLER
    [phiChain, sig2Chain, ~, rhoChain] = main_gibbs_rho(data, lambda, nSim);

    % PLOT SOME EVOLUTION GRAPH
    fileName = sprintf('Sig2Path_lambda_%4.2f', lambda);
    ExportSigma2Evolutiton(sig2Chain,fileName,rootPath);
    
    fileName = sprintf('RhoPath_lambda_%4.2f', lambda);
    ExportRhoEvolutiton(rhoChain,fileName,rootPath);

    % COMPUTE THE AVERAGE OF THE RESULTS
    [COV_EST, COVINV_EST] = estimationFromGibbsPaths(phiChain, sig2Chain, nBurnIn);
    [eLoss,qLoss] = getLoss(COV_EST, SIG_TRUE);

    % PLOT THE RESUILTS OF ESIMATION
    fileName = sprintf('Result_lambda_%4.2f.eps',lambda);
    ExportCovarianceMatricesGraph(SIG_TRUE, COV_EST, SIGINV, COVINV_EST, data, lambda, fileName, rootPath);

end

