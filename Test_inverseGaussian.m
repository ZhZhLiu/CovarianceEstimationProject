%% Test inverse gamma distribution
lambda = 1;
mu = 2;
size = 20000;
a = zeros(1,size);
for i = 1:size
    a(i) = ingrnd(mu,lambda);
end

% a = testFunction();