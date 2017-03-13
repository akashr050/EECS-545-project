function [ B ] = make_positive_definite(A)
%MAKE_POSITIVE_DEFINITE Summary of this function goes here
%   Detailed explanation goes here

i = 1;
[n,~] = size(A);
[~,p] = chol(A + i*speye(n));
while p > 0
    i= i+1;
    [~,p] = chol(A + i*speye(n));
end
A = A + i*speye(n);
B = A;
end

