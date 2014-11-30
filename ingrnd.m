function [ r ] = ingrnd(mu, lambda)

c = mu*randn(1)^2;
r = (mu/(2*lambda)) * (2*lambda + c - sqrt(4*lambda*c + c^2));
%invert = (rand(1).*(mu+r) > mu);
%r(invert) = mu(invert).^2 ./ r(invert);

test = rand(1)*(mu+r) > mu;
if(test)
    r = mu^2/r;
end

end

