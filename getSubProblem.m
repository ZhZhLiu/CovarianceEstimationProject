function [y,X] = getSubProblem(k, data)

y = data(:,k);
X = data(:,1:k-1);

end