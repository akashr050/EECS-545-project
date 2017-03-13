
% Author: Akash Rastogi 
% This code is to perform the analysis on ARA-Data
% The first section deals with the implementation of glasso algorithm and
% the subsequent section deals with PC algorithm 

%% Glasso Algorithm
clear;
% algo_identifier = 1;
% Generating the various kinds of synthetic data (cycle, chain, grid)
% Initialising the base parameters for data generation
% n = 100;
% p = 25;
% density = 0.80;
% [ X, H, G] = synthetic_data(n, p, @cycle_graph, density);
% 
rhos = 2.^(-(2:7)); 
ths  = [0.2 0.15 0.1 0.05 0.005];

load('AraData.mat');
[ ~, p] = size(X);
[ rho_best, th_best, ~] = cross_validation_ebic(H, X, rhos, ths, 5, 1);
JT_glasso_Gest = JT_UGMS(H, X, rho_best, th_best, 1);

H_glasso = ones(p, p);
[ rho_best, th_best, ~] = cross_validation_ebic(H, X, rhos, ths, 5, 2);
Glasso_Gest = get_estimate_glasso(H, X, rho_best, th_best);

subplot(2, 2, 1)
plot(graph(H));
title('H graph');

% Effect of thresholding on the covariance matrix
threshold_cov = (abs(cov(X)^(-1))>th_best).*(ones(p,p)-eye(p));
subplot(2, 2,2);
plot(graph(threshold_cov));
title('Effect of thresholding on Covariance matrix');

% G estimate from the JT-glasso algorithm
subplot(2, 2, 3)
plot(graph(JT_glasso_Gest));
title('G estimate from JT Glasso algorithm');

% G estimate from the Glasso algorithm
subplot(2, 2, 4)
plot(graph(Glasso_Gest));
title('G estimate from Glasso algorightm');

% save('aradata_analysis_cv_ebic');

%% PC algorithm
clear;

load('AraData.mat');
th = 0.2; % We used the threshold of 0.2 as per the vats paper

JT_PC_Gest = JT_UGMS(H, X, 0, th, 0);
PC_Gest = PC(X, 1, H, H, th) ;

subplot(2, 2, 1)
plot(graph(H));
title('H graph');

% Effect of thresholding on the covariance matrix
threshold_cov = (abs(cov(X)^(-1))>0.2).*(ones(p,p)-eye(p));
subplot(2, 2,2);
plot(graph(threshold_cov));
title('Effect of thresholding on Covariance matrix');

% G estimate from the JT-PC algorithm
subplot(2, 2, 3)
plot(graph(JT_PC_Gest));
title('G estimate from JT PC algorithm');

% G estimate from the PC algorithm
subplot(2, 2, 4)
plot(graph(PC_Gest));
title('G estimate from PC algorightm');

% save('aradata_analysis_pc');
