%Scripts to test the inverse gaussian distribution
mu = 1;
lambda = 2;


pd = makedist('InverseGaussian', 'mu', mu, 'lambda', lambda);
x = linspace(min(bins), max(bins), 1000);
y = pdf(pd,x);
hold on
plot(x,y);


size =10000;
results = invgaurnd(mu,lambda, size,1);

[counts, bins] = hist(results,100);
binWidth = bins(2) - bins(1);

plot(bins, counts/size/binWidth);

