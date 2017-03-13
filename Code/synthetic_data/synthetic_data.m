% Author: Yutong Wang

function [ X, H, G ] = synthetic_data(n, p, my_graph, density)
%SYNTHETIC_DATA Generate simulated data based on a graph model
%   Inputs: n = number of samples
%           p = number of features
%           my_graph = function handle of the base graph (note, the base
%           graph function may or may not be already PSD because we will
%           test it and make it so if it is not PSD already
%           density = percent of zeros in the error (density closer to 1
%           means less error, density closer to 0 means more error)
%   Outputs: X = n-by-p data matrix
%            H = the perturbed graph with noise
%            G = the true graph
    G = my_graph(p);
    [~, PSD] = chol(G);
    if PSD > 0 % if not PSD
        A = make_positive_definite(G);
    else % if is PSD
        A = G;
    end
    A = make_positive_definite(G);
    X = generate_X(A, zeros(1,p), n);
    noise = rand(p,p).*(ones(p,p)-eye(p));
    noise = (noise+noise')>(density*2);
    H = max((A>0)-eye(p), noise);
end

