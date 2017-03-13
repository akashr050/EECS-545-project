% Author: Yutong Wang 
% Adapted from: Xiaohui Chen (https://publish.illinois.edu/xiaohuichen/code/graphical-lasso/)

function [Theta, W] = glasso(S, rho, maxIt, tol)
% Solve the graphical Lasso
% minimize_{Theta > 0} tr(S*Theta) - logdet(Theta) + rho * ||Theta||_1
% Ref: Friedman et al. (2007) Sparse inverse covariance estimation with the
% graphical lasso. Biostatistics.
% Note: This function needs to call an algorithm that solves the Lasso
% problem. Here, we choose to use to the function *lassoShooting* (shooting
% algorithm) for this purpose. However, any Lasso algorithm in the
% penelized form will work.
%
% Input:
% S -- sample covariance matrix
% rho --  regularization parameter
% maxIt -- maximum number of iterations
% tol -- convergence tolerance level
%
% Output:
% Theta -- inverse covariance matrix estimate
% W -- regularized covariance matrix estimate, W = Theta^-1

    p = size(S,1);

    if nargin < 4, tol = 1e-6; end
    if nargin < 3, maxIt = 1e2; end

    % Initialization
    W = S + rho * eye(p);   % diagonal of W remains unchanged
    W_old = W;
    i = 0;

    % Graphical Lasso loop
    while i < maxIt,
        i = i+1;
        for j = p:-1:1,
            jminus = setdiff(1:p,j);
            W_11 = W(jminus,jminus);
            s_12 = S(jminus,j);
            b = lassoShooting(W_11, s_12, rho, maxIt, tol);
            W(jminus,j) = W_11 * b; % w_12
            W(j,jminus) = W(jminus,j)'; % w_12^T
        end
        % Stop criterion
        if norm(W-W_old,1) < tol, 
            break; 
        end
        W_old = W;
    end
    if i == maxIt,
        fprintf('%s\n', 'Maximum number of iteration reached, glasso may not converge.');
    end

    Theta = W^-1;
end

% Shooting algorithm for Lasso (unstandardized version)
function b = lassoShooting(V, u, lambda, maxIt, tol)
    % V = W11
    % u = s12
    % lambda = rho
    
    if nargin < 4, tol = 1e-6; end
    if nargin < 3, maxIt = 1e2; end

    % Initialization
    [p,~] = size(V);
    if p == 1 % dimension = 1 case
        b1 = (u+lambda)/V;
        b2 = (u-lambda)/V;
        if (1/2)*(sqrt(V)*b1 - sqrt(1/V)*u)^2 + lambda*sign(b1) < (1/2)*(sqrt(V)*b2 - sqrt(1/V)*u)^2 + lambda*sign(b2)
            b = b1;
        else
            b = b2;
        end
    else
        b= zeros(p,1);
        b_old = b;
        i = 0;
        ST = @(x,t) sign(x)*max(0,abs(x)-t); 
        % Shooting loop
        while i < maxIt,
            i = i+1;
            for j = 1:p,
        %       j = randi([1 p]);
                jminus = mysetdiff(1:p,j);
                b(j) = ST(u(j)-b(jminus)'*V(jminus,j),lambda)/V(j,j);
            end
            delta = norm(b-b_old,1);    % Norm change during successive iterations
            if delta < tol, break; end
            b_old = b;
        end
        if i == maxIt,
            fprintf('%s\n', 'Maximum number of iteration reached, shooting may not converge.');
        end
    end
end