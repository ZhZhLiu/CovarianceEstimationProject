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
    case 'diagDecay'
        sig2 = p:-1:1;
    case 'csym'
        for i = 2:p
          for j = i-1:-1:1
            T(i,j) = - rho / (1+ (i-1)*rho);
          end 
        end
        si = 1:p;
        sig2 = sig2 * (1- (si -1 ) * rho^2 ./ (1+(si-1)*rho));
end

TINVt =  T'\eye(p); % (T')^{-1}

% Generate the D matrix
DINV = diag(ones(1,p) ./ sig2);
D    = diag(ones(1,p) .* sig2);

% Based on the modified choleksy decomposition
SIGINV = T' * DINV * T;
SIG    = TINVt' * D * TINVt;
cholFact = sqrt(D) * TINVt;

% Simulation based on the cholesky factor
data = randn(n,p) * cholFact;
end

