

%% Parameter block
n = 100000;       % number of observations
p = 3;         % dimensions
sig2 = 0.01;     % variance of the data
rho = 0.8;      % exponential decay of the T matrix


cutoff = false;


% Generate the T matrix

T = eye(p);

if(~cutoff)

    for i = 2:p
       for j = i-1:-1:1
          T(i,j) = - rho^(i-j);
       end 
    end

else
    for i = 2:p
        T(i,i-1) = -rho;
    end
end

TINVt = T'\eye(p);

DINV = eye(p) / sig2;
D = eye * sig2;


%% Based on the modified cholesky decomposition

SIGINV = T' * DINV * T;
SIG = TINVt' * D * TINVt;
cholFact = sqrt(sig2) * eye(p) * TINVt;

%% Simulation based on  the results

data = randn(n,p) * cholFact;