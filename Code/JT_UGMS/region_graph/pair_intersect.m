% Author: Vyas Ramasubramani

function [intersections, iter1, iter2] = pair_intersect(sets)
% Function to compute pairwise set membership for elements of sets
% Input: sets - A binary matrix where each row is a set and each column is
%               a possible element
% Output: intersections - A binary matrix where for each pair of rows 
%                           (i,j) in sets the row (i-1)*n + j in
%                           intersections contains the bitwise intersection
%                           of the rows
%         indices - A matrix whose (i, j) entry indicates the location of
%                   the intersection of rows i and j from sets in the
%                   output matrix intersections
% Notes: Does not include empty intersections, so using the iter1 and iter2
% outputs is necessary

n = size(sets, 1);
iter1 = repelem(1:n, n);
iter2 = repmat(1:n, 1, n);

% Remove nonintersecting pairs to minimize computation/memory use
combinations = (sets*sets')' > 0;
iter1 = iter1(combinations(:) > 0);
iter2 = iter2(combinations(:) > 0);

% Avoid duplicates due to symmetry (i, j) vs (j, i)
indices = iter2 > iter1;
iter1 = iter1(indices)';
iter2 = iter2(indices)';

check1 = sets(iter1, :);
check2 = sets(iter2, :);

% Element-wise multiply tells us membership since these are bits
intersections = check1.*check2;


% 
% n = size(sets, 1);
% iter1 = repelem(1:n, n);
% iter2 = repmat(1:n, 1, n);
% indices = iter2 > iter1;
% iter1 = iter1(indices)';
% iter2 = iter2(indices)';
% 
% check1 = sets(iter1, :);
% check2 = sets(iter2, :);
% 
% % Element-wise multiply tells us membership since these are bits
% intersections = check1.*check2;
end