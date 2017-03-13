% author: Akash Rastogi, Yutong Wang
function [ rho_best, th_best, loss_matrix ] = cross_validation_ebic(H, X, rhos, ths, B, algo_id)
%CROSS_VALIDATION Summary of this function goes here
%   input:      rhos        - a row of rhos to test
%               ths         - a row of thresholds to test
%               B    - number of bins you want
%               algo_id - is the identifier for algorithm we want to investigate
%                         JT_glasso: 1
%                         Glasso: 2
%   outputs:    rho         - best rho
%               threshold   - best threshold
%               loss_matrix - matrix of loss for the given rhos and ths
%               combinations
    [n, p]= size(X);
    L = floor(n/B);
    bins = zeros(B, n);
    
    for b = 1:B-1
        bins(b,(b-1)*L+1:b*L) = 1;
    end
    bins(B,(B-1)*L+1:end) = 1;
    bins = ones(B,n) - bins;
    bins = bins == 1;
    rho_best = rhos(1);
    th_best = ths(1);
    loss_best = Inf;
    R = numel(rhos);
    T = numel(ths);
    loss_matrix = zeros(R,T);
    for r = 1:R
    for t = 1:T
        rho = rhos(r);
        th  = ths(t);
        trivial = 1; % unless proven otherwise, a certain combination of rho and threshold produces all zeros
        loss = 0;
        for b1 = 1:B
            [obs, p] = size(X(bins(b1,:),:));
            if algo_id == 1
                Gest = JT_UGMS(H, X(bins(b1,:),:), rho, th, algo_id);
            elseif algo_id == 2
                Gest = get_estimate_glasso(H, X(bins(b1,:),:), rho, th);
            end
            if sum(Gest(:)) > 0 % if any Gest1 nontrivial, then it's not all zeros
                trivial = 0;
            end
%           Implementation of EBIC
            e_gest = sum(Gest(:)/2);
            Gest = Gest + eye(p);
            EBIC = obs * (log(det(Gest)) - trace(cov(X) * Gest)) + e_gest * log(obs) + 4*0.5* e_gest * log(p);                
            EBIC = - EBIC;
            loss = loss + EBIC;
        end
        loss_matrix(r,t) = loss;
        % if not trivial and loss improved...
        if (trivial==0) && (loss < loss_best)
            loss_best = loss;
            th_best = th;
            rho_best = rho;
        end
    end
    end
    if loss_best == Inf
        error('Cross Validation failed, all estimated graphs were trivial');
    end
end

