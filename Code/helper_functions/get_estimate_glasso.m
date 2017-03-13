function [iGestimate] = get_estimate_glasso(HR, X, rho, threshold)
    [Theta, ~] = glasso(cov(X), rho);

    %     disp('hr');
%     disp(hr);
%     disp('iGestimate');
%     disp(iGestimate);
    iGestimate = Theta .* HR;
    iGestimate = abs(iGestimate) > threshold;
    iGestimate = triu(iGestimate,1)+tril(iGestimate,-1);
end
    