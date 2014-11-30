%% Regression testing of the density for pdf
clear
testN = 100;
testResult = zeros(testN,2);

mus = abs(randn(testN,1))*2;
lams = abs(randn(testN,1))*2;
xs  = abs(randn(testN,2));


k = 10;

sig2s = abs(randn(testN,1));
phis = randn(testN,k-1);
tau2invs = abs(randn(testN,k-1));
rhos = rand(testN,1);

j=6;

for g=1:10000

for i = 1:testN
    
    x = xs(i,:);
    phi = phis(i,:);
    tau2inv = tau2invs(i,:);
    sig2 = sig2s(i,:);
    lamda = 2;
    rho = rhos(i);
    
    p1_1 = tauinv2Density(x(1),j, phi, tau2inv, sig2, rho, lamda, k);
    p1_2 = tauinv2Density(x(2),j, phi, tau2inv, sig2, rho, lamda, k);
    
    p2_1 = tauinv2Density2(x(1), j, k, rho, phi, tau2inv, sig2,lamda);
    p2_2 = tauinv2Density2(x(2), j, k, rho, phi, tau2inv, sig2,lamda);
    
    testResult(i,1) = p1_1/p1_2;
    testResult(i,2) = p2_1/p2_2;

end

end