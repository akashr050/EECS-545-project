% Author: Yutong Wang, Akash Rastogi, Yiqian Gan

function [rBar] = get_RBar(i, L, regions, edges)
% Code to get Rbar
% Input: index - which clique
%           L - which layer
% Output: rBar - is union of ancestors values as per the vats paper   

    children = i; % initialize to just the current region
    for layer = L:-1:2
        parents = edges{layer-1}(:,1);
        children = parents(ismember(edges{layer-1}(:,2),children));
    end
    rBar = find(max(regions{1}(children,:),[], 1));
end
