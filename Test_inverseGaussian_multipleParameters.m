%% TEST MULTIPLE PARAMETERS GENERATION
clear;
mus = 1:10;
nObs = 10000;
lams = mus*2;
r = invGaussianMultipleParams(lams,mus, nObs);

%%

for i = 1:10
    [counts, bins] = hist(r(:,i),100);
    binWidth = bins(2) - bins(1);
    plot(bins, counts/nObs/binWidth);
    hold on;
    
    pd = makedist('InverseGaussian', 'mu', mus(i), 'lambda', lams(i));
    x = linspace(min(bins), max(bins), 1000);
    y = pdf(pd,x);
    hold on
    plot(x,y);
    
end

