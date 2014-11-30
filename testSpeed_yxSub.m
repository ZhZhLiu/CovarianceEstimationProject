for i = 1:10000
k = 20;
[y,X] = getSubProblem(k,data);

Subproblem{k}.y;
Subproblem{k}.X;
end