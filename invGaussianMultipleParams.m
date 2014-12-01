function [ r ] = invGaussianMultipleParams(lams,mus, nObs)


nParams = length(mus);
muMat   = ones(nObs, nParams) * diag(mus);
lamMat  = ones(nObs, nParams) * diag(lams);

c = randn(nObs, nParams).^2 .* muMat;
r = (muMat ./ (2*lamMat)) .* (2 * lamMat + c - sqrt( 4* lamMat .* c + c.^2));
test = rand(nObs, nParams) .* (muMat + r) > muMat;
r = (1-test).* r + test .* ( muMat.^2 ./ r);

end

