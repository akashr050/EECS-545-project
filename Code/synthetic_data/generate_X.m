function [ X ] = generate_X(A,mu,n)
%GENERATE_X make synthetic data
% A is the weighted adjacency matrix (without 
% requires mu to be a 1-by-p matrix
% n is the number of samples
% X is n by size(mu)
    Sigma = inv(make_positive_definite(A));
    X= mvnrnd(repmat(mu,n,1),Sigma); % data amatrix of size n-by-p
end

