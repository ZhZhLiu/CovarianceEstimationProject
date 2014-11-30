function [data,SIGINV, SIG, T, DINV] = generateDataset(n,p,sig2,rho,type)
% This function is used for generating the dataset for simulation

% Generate the T matrix
T = eye(p);
switch type
    case 'cutoff'
        for i = 2:p
            T(i,i-1) = -rho;
        end
    case 'decay'
        for i = 2:p
          for j = i-1:-1:1
            T(i,j) = - rho^(i-j);
          end 
        end
    case 'step'
        for i = 2:p
          for j = i-1:-1:1
            T(i,j) = max(0, 1- rho*(i-j));
          end
        end
end

TINVt =  T'\eye(p); % (T')^{-1}

% Generate the D matrix
DINV = eye(p) / sig2;
D    = eye(p) * sig2;

% Based on the modified choleksy decomposition
SIGINV = T' * DINV * T;
SIG    = TINVt' * D * TINVt;
cholFact = sqrt(sig2) * eye(p) * TINVt;

% Simulation based on the cholesky factor
data = randn(n,p) * cholFact;
end

