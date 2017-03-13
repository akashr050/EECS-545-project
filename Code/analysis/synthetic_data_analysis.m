% Author: Akash Rastogi 
% This code was used to perform the ad-hoc analysis on the synthetic data

%% Glasso analysis 
clear;
algo_identifier = 1;
% Generating the various kinds of synthetic data (cycle, chain, grid)
% Initialising the base parameters for data generation
n = 100;
p = 25;
density = 0.75;
[ X, H, G] = synthetic_data(n, p, @cycle_graph, density);

rhos = 2.^(-(2:12)); 
ths  = [0.2 0.15 0.1 0.05 0.005 0.0005];
metric = zeros(length(rhos) * length(ths), 6);

counter = 1;
for i = 1:length(rhos)
    for j = 1:length(ths)
        % For JT_glasso alorightm uncomment the line below and comment out
        % the three lines below that
        % Gest = JT_UGMS(H, X, rhos(i), ths(j), algo_identifier);
        
        [ ~, p] = size(X);
        H = ones(p, p);
        Gest = get_estimate_glasso(H, X, rhos(i), ths(j));
        met = metrics(Gest, G); 
        if (met(2) ~= 0)
            my_metric = (1 - met(1)) * met(2);
        else
            my_metric = met(1);
        end
        metric(counter, :) = [met my_metric rhos(i) ths(j)];
        counter = counter+1;
    end 
end

[~, best_index] = max(metric(:,4));
best_parameters = metric(best_index,:);

%% PC analysis


clear;
% Generating the various kinds of synthetic data (cycle, chain, grid)
% Initialising the base parameters for data generation
n = 100;
p = 25;
density = 0.75;
[ X, H, G] = synthetic_data(n, p, @cycle_graph, density);

% Gest = JT_UGMS(H, X, 0, 0.2, 0);
Gest = PC(X, 1, H, H, 0.2);
met = metrics(Gest, G); 
met = [met (1-met(1)) * met(2)];