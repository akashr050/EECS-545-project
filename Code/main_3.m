% Author: Yiqian Gan, Vyas Ramasubramani, Akash Rastogi, Yutong Wang
% This file shows the implementation of our JT-glasso algorithm.
% Just plug in the dataset and play 

% Loading the dataset
load('AraData.mat');

% Initializing the parameters 
rhos = 2.^(-(2:7)); 
ths  = [0.2 0.15 0.1 0.05 0.005];
CV_bins = 5;
algo_identifier = 1; % This is an identifier for glasso algorithm

% Using cross validation to tune the parameter: rho and threshold 
[ rho_best, th_best, loss_matrix] = cross_validation_ebic(H, X, rhos, ths, CV_bins, 1);

% Using the optimised parameters to estimate graph
JT_glasso_Gest = JT_UGMS(H, X, rho_best, th_best, 1);


%% plotting

% H graph
subplot(1, 2, 1)
plot(graph(H));
title('H graph');

% G estimate from JT-Glasso algorithm
subplot(1, 2, 2)
plot(graph(JT_glasso_Gest));
title('G estimate from JT Glasso algorithm');
