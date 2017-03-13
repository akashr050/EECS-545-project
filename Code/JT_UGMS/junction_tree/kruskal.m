% Author: Vyas Ramasubramani

function MST = kruskal(A)
%kruskal Apply Kruskal's algorithm to find the minimal spanning tree of A,
%an undirected graph

% Possible errors:
% A is not 2d; A is not square; A is not symmetric
% A must be an undirected graph
if numel(size(A)) > 2
    error('The input must be a 2d array')
elseif size(A, 1) ~= size(A, 2)
    error('The input must be a square matrix')
elseif ~isequal(A, A')
    error('The matrix input must be symmetric! This function does not support undirected graphs')
end

weight_matrix = tril(A, -1); % Symmetry lets us only deal with half the matrix
weight_matrix(weight_matrix == 0) = Inf; % Make min ignore zeros
tree_vector = 1:size(A, 1); % Vector indicating which tree each vertex belongs to
edge_matrix = zeros(size(A)); % Binary matrix indicating when two elements are connected

% The termination conditions are either that we have constructed a spanning
% tree, or we have exhausted all edges
% If we run out of edges as the terminating condition (rather than putting
% everything into one tree), what does that mean? Is that even a viable
% scenario?
while numel(unique(tree_vector)) ~= 1 && ~all(all(weight_matrix == Inf))
    % Pick the smallest edge. Ties are by broken by choosing the smallest
    % index, which is the min function's default behavior
    [M, I] = min(weight_matrix);
    [~, col_index] = min(M);
    row_index = I(col_index);

    % Eliminate this edge from the weight matrix
    weight_matrix(row_index, col_index) = Inf;
    
    % If new edge would connect two previously disjoint trees, add it and
    % merge the trees
    if tree_vector(row_index) ~= tree_vector(col_index)
        edge_matrix(row_index, col_index) = 1;
        
        % For consistency, always move towards lower indices
        new_index = min(tree_vector([row_index, col_index]));
        tree_vector(tree_vector == tree_vector(row_index)) = new_index;
        tree_vector(tree_vector == tree_vector(col_index)) = new_index;
    end
end

% We could equivalently just return the upper (or lower) triangular part.
MST = (edge_matrix + edge_matrix').*A;
end