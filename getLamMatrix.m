function [lamInv] = getLamMatrix(rho, dim)

lamInv = eye(dim) * (1+rho^2);
lamInv(1,1) = 1;
lamInv(dim,dim) = 1;

offDiagMat = zeros(dim);

for i = 1:dim-1
    offDiagMat(i,i+1) = -rho;
end
lamInv = lamInv + offDiagMat + offDiagMat';

lamInv = lamInv / (1-rho^2);
end

