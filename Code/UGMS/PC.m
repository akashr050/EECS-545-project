% Author: ??? Yiqian???

function G = PC(X,k,H,L, xi)
% Implements the PC algorithm
% Inputs:
% X: n i.i.d. observations. X_rbar
% k: An integer that controls the computational complexity of PC.
% H: A graph that contains all the true edges G* ->HUG_rbar
% L: A graph that contains the edges that need to be estimated. ->H_r'
% Output: A graph Ghat that contains edges in L that are estimated to be in G?

% Initialization
p = size(X,2); % number of nodes;
n = size(X,1);
% xi = 0.2; % 0.2 thresholding parameter
Delta = 0:k; 
pcType = 1;
xiChosen = 0;
displayInd = 0;


if isempty(H)
    H = ones(p,p);
end

if isempty(L)
    L = H;
end
% Initialization
G = L;

if min(Delta) == 0
    cc = mycorr(X);
    bool = (abs(cc) < xi);
    G(bool == 1) = 0;
    H(bool == 1) = 0;
    Delta = setdiff(Delta,0);
end

for k = Delta
    % find all edges that need to estimated
    ind_edges = find(triu(G) == 1);
    %ind_edges = ind_edges(randperm(length(ind_edges)));
    num_edges = length(ind_edges);
    for e = 1:num_edges
        [i j] = ind2sub([p p], ind_edges(e)); % the edge
%         nei = FindNeighborsUndirected(H,i,1);
        % find the nearest neighbor
        nei_my_list = H(i,:)==1;
        nei_my_list(i)=0;
        temp = 1:p;
        nei = temp(nei_my_list);
        
        % find the nearest neighbor
        nej_my_list = H(j,:)==1;
        nej_my_list(j)=0;
        temp = 1:p;
        nej = temp(nej_my_list);

        Sij = SeparatorSearchSpaceNew(nei,nej,i,j,H,pcType,k);
%         Sij = SeparatorSearchSpace(nei_my,nej_my,i,j,H,pcType,k);
        if ~isempty(Sij)
            comb = combnk(Sij,k);
            for ct = 1:size(comb,1)
                S = comb(ct,:);
                % conditional independence
                % Can replace this by PartialCorrCoef(cov(X),i,j,S,n) too
                cc = partialcorr(X(:,i),X(:,j),X(:,S));
                bool = (abs(cc) < xi); %bool = 1 => X_i indep X_j | X_S
                if (bool == 1)
                    G(i,j) = 0; G(j,i) = 0; H(i,j) = 0; H(j,i) = 0;
                    break;
                end
            end
        else
            H(i,j) = 0; H(j,i) = 0;
        end
    end
end

% G = setdiag(G,0);
G = G - diag(diag(G));
end