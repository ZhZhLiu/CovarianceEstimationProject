%% Regression testing of the density for pdf
testN = 100;
testResult = zeros(testN,2);

mus = abs(randn(testN,1))*2;
lams = abs(randn(testN,1))*2;

xs  = abs(randn(testN,1));

for i = 1:testN
    mu = mus(i);
    lam = lams(i);
    x = xs(i);
    
    pd = makedist('InverseGaussian','mu',mu,'lambda', lam);
    testResult(i,1) = pdf(pd,x);
    testResult(i,2) = invGauPdf(x,lam,mu);
   
end