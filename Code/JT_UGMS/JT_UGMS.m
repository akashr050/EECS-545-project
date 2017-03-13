% Author: Akash Rastogi    

function [Gest] = JT_UGMS(H, X, rho, threshold, algo)
% Implements the overall Junction Tree UGMS framework
% Algo: Identifier for algorithm implemented
%        1: Glasso; 0 - PC

    elim_order = get_elimination(H);
    H_triangulated = fill_graph(H, elim_order);
    cliques = gen_maximal_cliques(H_triangulated);
    JG = get_JG_from_mcliques(cliques);
    JT = -kruskal(-JG);
    [regions, edges] = region_graph(JT, cliques);
    % for i = 1:length(RegionMap)
    %     disp(i)
    %     disp(RegionMap(i))
    % end

    [n, ~] = size(H);
    Gest = zeros(n, n);
    EET = zeros(n,n);
    HUG = H;

    while(1)
        L = 1; % layer number
        flag = 1; % flag = 0 if the layer has a non-empty HR
        EET = zeros(n, n);
        while(flag) % determine which layer you're running on, if your layer has a non-empty HR, then break out of this while loop
            [m, ~] = size(regions{L});
            for i = 1:m % loops over regions in the layer
                if L ==length(regions) % on the last layer, pass in no children because they all have been ran over
                    [HR, EET] = HR_graph(i, H, regions{L}, zeros(0,n), zeros(0,2), EET); % edge case for the bottom layer
                else
                    [HR, EET] = HR_graph(i, H, regions{L}, regions{L+1}, edges{L}, EET); % non-edge case
                end

                if(sum(HR(:)) > 0)  % if all HR graphs are empty (i.e. no edges), then move down a layer
                    RBar = get_RBar(i, L, regions, edges);   % self + parents
                    
                    iX = X(:, RBar);
                    iHUG = HUG(RBar, RBar);
                    if algo == 1     % glasso
                        iGest = get_estimate_glasso(HR(RBar,RBar), iX, rho, threshold);
                    elseif algo == 0 % PC
                        iGest = get_estimate_PC(iX, 1, iHUG, HR(RBar,RBar), threshold);
                    end
                    Gest(RBar, RBar) = Gest(RBar, RBar) + iGest; 
 
                    
                    
%                     % For Glasso, uncomment the line below 
%                      % iGest = get_estimate_glasso(HR(RBar,RBar), iX, rho, threshold);
%                     % For PCA run the code below
%                      iHUG = HUG(RBar, RBar);
%                       iGest = get_estimate_PC(iX, 1, iHUG, HR(RBar,RBar), threshold);
%                     Gest(RBar, RBar) = Gest(RBar, RBar) + iGest; 
                    
                    flag = 0; % if any HR is non-empty, set flag to 0
                end
            end
            H = max(H - EET, 0);
            if flag == 1
                if L == length(regions) % if you hit the bottom layer
                    break
                end
                L = L+1; % go down a layer
            end    
        end
        if sum(H(:)) == 0 % make sure that H is non-empty, otherwise, we're done
            break
        end
        
        % reconstruct the region graph
        HUG = H + Gest;
        elim_order = get_elimination(HUG);
        A_triangulated = fill_graph(HUG, elim_order);
        cliques = gen_maximal_cliques(A_triangulated);
        JG = get_JG_from_mcliques(cliques);
        JT = -kruskal(-JG);
        [regions, edges] = region_graph(JT, cliques);
    end
end

