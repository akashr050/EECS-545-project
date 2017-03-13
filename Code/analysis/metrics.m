% Author: Akash Rastogi
function [metrics] = metrics(g, gStar)
% Calculate the performance metrics for the estimated graph
% Metrics is a vector of fdr, tpr and ed. Description of all these terms is
% provided in the report

    g = triu(g);
    gStar = triu(gStar);
    fdr = sum(sum(max(g - gStar, 0))) / sum(g(:));
    tpr = sum((g(:) == gStar(:))' * g(:)) / sum(gStar(:));
    ed = sum(sum(max(abs(g - gStar), 0)));
    metrics = [fdr tpr ed];
end