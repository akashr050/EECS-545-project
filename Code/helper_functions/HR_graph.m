% Author: Yutong Wang, Akash Rastogi, Yiqian Gan

function [HR, EET] = HR_graph(i, H, parent_layer, children_layer, edges, EET)
% Code to get Hr
% Input: RegionMap - is the cluster map of region graph
%        RG - is the region graph
%        index - is the index corresponding to which we require hr graph
%        A - is the original H graph from which we will derive hr
%           EET - is the n-by-by adjacency matrix of estimated edges
% Output: hr - is the hr graph as per its defination by vats paper 

% IMPORTANT NOTE: for the last layer, call hrGraph3(i, H, parent_layer, zeros(0,n), zeros(0,2), EET)
% because obviously there is no more children, they all have been ran over

    [n,~] = size(H);
    HR = zeros(n,n);
    R = parent_layer(i,:)==1; % region R^l_i logical vector
    childrens = edges(:,2);
    childrens_of_R = childrens(edges(:,1)==i); % find children of R
    HR(R,R) = H(R,R);
    for c =childrens_of_R' % matlab loops over the rows only!
        child = children_layer(c,:)==1;
        HR(child,child) = 0;
    end
    EET = EET + HR;
end
