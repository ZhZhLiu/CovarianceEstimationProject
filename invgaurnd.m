function [r] = invgaurnd(mu, lambda, m,n)

c = mu*randn(m,n).^2;
r = (mu/(2*lambda)) * (2*lambda + c - sqrt(4*lambda.*c + c.^2));

test = rand(m,n).*(mu+r) > mu;

r = (1-test).* r + test .* ( mu^2 ./ r);

end

