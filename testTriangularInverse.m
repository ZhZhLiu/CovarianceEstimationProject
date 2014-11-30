A = randn(200);
U = triu(A);
I = eye(200);
inv(A);
inv(U);
U\I;