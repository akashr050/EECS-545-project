% Author: Vyas Ramasubramani

function [regions, edges] = region_graph(JT, clusters)
% Function to generate the Region map from the given junction tree
% Input: JT - Adjacency matrix for junction tree
%        clusters - Clusters that are vertices in junction tree
% Output: regions - Cell list where each element represents a layer. Each
%                   element is a matrix whose rows are the regions in that 
%                   layer, and the columns are binary indicating which
%                   elements of your original graph are in that region
%         edges - Cell list where each element represents edges between one
%                 layer and the next. Each element is an nx2 matrix that
%                 just contains the source and sink of the edge

% NOTE: CURRENTLY THE BOTTLENECK OPERATION IS THE UNIQUE, WHICH INVOLVES
% SORTING ROWS. I HAVEN'T BOTHERED TO TRY AND FIND A FAST WAY TO DO IT
% SINCE I DOUBT IT WILL BE A REAL PROBLEM, BUT IF IT IS THAT'S WHAT WE NEED
% TO OPTIMIZE. WE COULD EASILY DO SOMETHING BY JUST PAIRWISE MULTIPLYING
% SIMILAR TO MY PAIR_INTERSECT FUNCTION, BUT THE TRICK WILL BE FINDING A
% WAY TO ALSO GET BACK INDICES SINCE I USE THAT PART OF THE UNIQUE OUTPUT
% AS WELL RIGHT NOW

% Since we don't know the depth we can just use cell lists. Might be good
% to preallocate and just double in size every time it's not big enough
regions = {};
edges = {}; % List of nx2 matrices. First column is region in current level, second column is region in next level

% Special logic for the first two levels
regions{1} = clusters;
n = size(JT, 1);
tocheck = triu(JT);
indices = [tocheck(:), repelem(1:n, n)', repmat(1:n, 1, n)'];
ii = indices(indices(:, 1) > 0, 2);
ij = indices(indices(:, 1) > 0, 3);
check1 = regions{1}(ii, :);
check2 = regions{1}(ij, :);
intersections = check1.*check2;
[regions{2},~,ic] = unique(intersections,'rows');
edges{1} = [ii, ic;ij, ic];

% We might already be done.
if all(sum(regions{1}, 2) <= 2)
    return
end

% Now we run a while loop until the termination condition is reached
i = 2;
while 1
    % I fixed memory issues for the initial case, but it still happens
    % later just because eventually the number of possible intersections is
    % too large. If it's necessary we might have to partition the
    % intersection logic to intersect a subset, remove empty intersections,
    % then keep doing that in a loop 
    [intersections, ii, ij] = pair_intersect(regions{i});
    check_intersect = sum(intersections, 2) > 1; % Avoids cardinality 1 intersections
    ii = ii(check_intersect, :);
    ij = ij(check_intersect, :);
    intersections = intersections(check_intersect, :);

    [regions{i+1},~,ic] = unique(intersections,'rows');
    edges{i} = [ii, ic;ij, ic];
    
    % Termination condition
    if all(sum(regions{i+1}, 2) <= 2)
        % In case the last layer has all regions of size exactly 2, then we
        % will have added an empty layer so remove it if needed
        if all(sum(regions{i+1}, 2) == 0)
            regions = regions(~cellfun('isempty', regions));
            edges = edges(~cellfun('isempty', edges));
        end
        break
    end

    i = i + 1;
end
end