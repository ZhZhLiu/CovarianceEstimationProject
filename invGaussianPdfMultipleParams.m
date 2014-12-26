function [pdf] = invGaussianPdfMultipleParams(obs, mus, lams)

[nObs, nParams] = size(obs);
muMat   = ones(nObs, nParams) * diag(mus);
lamMat  = ones(nObs, nParams) * diag(lams);

part1 = (lamMat ./ ( 2* pi * obs.^3) ).^(0.5);
part2 = exp( - lamMat.* (obs -muMat).^2 ./ ( 2 .* muMat.^2 .* obs));

pdf = part1 .* part2;

end

