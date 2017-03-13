function [A,i] = generatesparseSPDmatrix(n,density)
% Generate a sparse n x n symmetric, positive definite matrix with
%   approximately density*n*n non zeros

A = abs(sign(sprandsym(n,density))); % generate a random n x n matrix

% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI
i = 1;
[~,p] = chol(A + i*speye(n));
while p > 0
    i= i+1;
    [~,p] = chol(A + i*speye(n));
end
A = A + i*speye(n);
end
