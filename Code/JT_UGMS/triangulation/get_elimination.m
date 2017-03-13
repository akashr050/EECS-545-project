%Author: Vyas Ramasubramani

function ordering = get_elimination(A)
%get_elimination Computes a perfect elimination ordering on A using the
%greedyDegree heuristic
% Inputs:
%   A: A graph
%
% Output: A vector containing an elimination ordering on A

% Possible errors:
% A is not 2d; A is not square; A is not symmetric
% A must be an undirected graph
if numel(size(A)) > 2
    error('The input must be a 2d array')
elseif size(A, 1) ~= size(A, 2)
    error('The input must be a square matrix')
elseif ~isequal(A, A')
    error('The matrix input must be symmetric! This function does not support directed graphs')
end

n = size(A, 1);
H = tril(A, -1) + triu(A, 1); % Ignore diagonal
H = H > 0; % We can't do anything with weights for this algorithm, so just get rid of them
ordering = zeros(1, n);
checked = zeros(1, n);

for i = 1:n-1
    degrees = sum(H);
    % To handle initially disconnected graphs
    [~, v] = min(degrees+checked); % By default, breaks ties by earlier index
    
    % Determine which vertices to connect to make a clique and join them
    vertices = find(H(:, v));
    L = length(vertices);
    H(vertices,vertices) = ones(L,L) - eye(L);
    
    H(:, v) = 0;
    H(v, :) = 0;
    ordering(i) = v;
    checked(v) = Inf;
end

ordering(n) = (n*(n+1))/2 - sum(ordering);
end