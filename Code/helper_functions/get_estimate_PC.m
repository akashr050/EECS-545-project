function [iGestimate] = get_estimate_PC(X, k, HUG, HR, threshold)
    iGestimate = PC(X, k, HUG, HR, threshold) .* HR;
end

