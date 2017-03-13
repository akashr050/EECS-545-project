%Author: Vyas Ramasubramani

function A_triangulated = fill_graph(A, elim_order)
%fill Uses an elimination ordering to triangulate A
% Inputs:
%   A: A graph
%   elim_order: An elimination ordering on A
%
% Output: The triangulated version of A

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

n = length(elim_order);
A_triangulated = A;
for i = 1:n
    v = elim_order(i);
    remaining_vertices = elim_order(i:n);
    neighbors = A_triangulated(v, remaining_vertices).*remaining_vertices;
    neighbors = neighbors(neighbors ~= 0);
    s = length(neighbors);
    A_triangulated(neighbors, neighbors) = ones(s, s) - eye(s);
end

end