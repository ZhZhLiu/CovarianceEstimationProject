function [entropyLoss, quadLoss] = getLoss(EST_COV, TRUE_COV)
% The function is for measuring the loss of the estimator with resepect to
% the true covariance matrix..

RATIO = EST_COV \ TRUE_COV;
n = length(TRUE_COV);

entropyLoss = trace(RATIO) - log(det(RATIO)) - n;
quadLoss = (trace(RATIO - eye(n))).^2;



end
