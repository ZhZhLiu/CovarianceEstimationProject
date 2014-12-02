function [COV_EST, COVINV_EST] = estimationFromGibbsPaths(phiChain, sig2Chain, nBurnIn)

    [p,~] = size(sig2Chain);

    T_EST = -mean(phiChain(:,:,nBurnIn:end),3) + eye(p);
    S_EST = mean(sig2Chain(:,nBurnIn:end), 2);
    DINV_EST = diag(1./S_EST);

    COVINV_EST = T_EST' * DINV_EST * T_EST;
    COV_EST = inv(COVINV_EST);
end

