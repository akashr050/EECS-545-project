% Author: Vyas Ramasubramani

function cliques = gen_maximal_cliques(A)
%gen_maximal_cliques apply the Bron-Kerbosch algorithm to find maximal
%cliques on A
% Inputs:
%   A: A graph
%
% Output: A matrix whose rows are binary vectors indicating the clique
% membership of the feature in that column

% Generate degeneracy ordering on A
degrees = sum(A);
degree_matrix = sortrows([-degrees', (1:length(degrees))']);
deg_order = degree_matrix(:, 2)';

% Setup variables
num_features = size(A, 1);
P = 1:num_features; 
R = [];
X = [];
cliques = zeros(1, num_features);
num_cliques = 0;

for v = deg_order
    % Find neighbors of v
    neighbors_v = find(A(v, :));
    Bron_Kerbosch(myunion(R, v), myintersect(P, neighbors_v), myintersect(X, neighbors_v));
end
    
    function [] = Bron_Kerbosch(R, P, X)
        if(isempty(P) && isempty(X))
            num_cliques = num_cliques + 1;
            clique = zeros(1, num_features);
            clique(R) = 1;
            cliques(num_cliques, :) = clique;
        end
        
        % Use Koch's pivot selection criterion: the pivot u is the element
        % in the union of P and X such that the cardinality of the union of
        % P and the neighbors of u is maximized
        possible_pivots = myunion(P,X);
        in_P = zeros(num_features, 1);
        in_P(P) = 1;
        pcounts = A(possible_pivots,:)*in_P;
        [~,index] = max(pcounts);
        u = possible_pivots(index);
        
        % Now recurse
        neighbors_u = find(A(u, :));
        for w = mysetdiff(P, neighbors_u)
            neighbors_w = find(A(w, :));
            Bron_Kerbosch(myunion(R, w), myintersect(P, neighbors_w), myintersect(X, neighbors_w));
            P = mysetdiff(P, w);
            X = myunion(X, w);
        end
    end

cliques = unique(cliques, 'rows');
num_cliques = size(cliques, 1);
end