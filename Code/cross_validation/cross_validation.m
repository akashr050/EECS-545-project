% author: Yutong Wang
function [ rho_best, th_best, loss_matrix ] = cross_validation(H, X, rhos, ths, B, algo_identifier)
%CROSS_VALIDATION Summary of this function goes here
%   input:      rhos        - a row of rhos to test
%               ths         - a row of thresholds to test
%               B    - number of bins you want
%               my_UGMS - is the function handle for the algo we are
%               investigating
%   outputs:    rho         - best rho
%               threshold   - best threshold
    [n,~]= size(X);
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
    tic;
    h = sum(H(:));
    for r = 1:R
    for t = 1:T
        rho = rhos(r);
        th  = ths(t);
        trivial = 1; % unless proven otherwise, a certain combination of rho and threshold produces all zeros
        loss = 0;
        for b1 = 2:B
            Gest1 = JT_UGMS(H, X(bins(b1,:),:), rho, th, algo_identifier);
            if sum(Gest1(:)) ~= 0 || sum(Gest1(:)) ~= h % if any Gest1 nontrivial, then it's not all zeros
                trivial = 0;
            end
            for b2 = 1:b1
                Gest2 = JT_UGMS(H, X(bins(b2,:),:), rho, th, algo_identifier);
% Implementation of EBIC
% e_gest = sum(Gest(:)/2);
% [n, p] = size(X);
% Gest = Gest + eye(p);
% log(det(Gest)) - trace(cov(X) * Gest) + e_gest * log(n) + 4*0.5* e_gest * log(p)                
                loss = loss + sum(Gest1(:)~=Gest2(:))/(sum(Gest1(:)+Gest2(:))+1);
            end
        end
        rho
        th
        loss
        toc
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

