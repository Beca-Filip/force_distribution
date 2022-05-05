close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Create model
[model, vars] = form_casadi_model_normalized();
model.solver('ipopt');

% Create weight vector
alpha = zeros(17, 1);
alpha(2) = 1;

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial_list = 9;
speed_list = 5;
leg_list = 1;
sample_list = 1:101;

% Calculate Forces
Fout = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

%% Compare

figure;
hold all;
plot_vector_quantities(sample_list-1, squeeze(data.f(:,:,sample_list,trial_list,speed_list,leg_list)), [], 'LineWidth', 2, 'Color', [0,0,0]);
plot_vector_quantities(sample_list-1, squeeze(Fout), [], 'LineWidth', 2, 'Color', [.8,0,0]);
sgtitle('Data vs. Minimum Norm solution.');

figure;
hold all;
plot(sample_list-1, sum(squeeze(data.f(:,:,sample_list,trial_list,speed_list,leg_list)).^2, 1), 'LineWidth', 2, 'Color', [0,0,0]);
plot(sample_list-1, sum(squeeze(Fout).^2, 1), 'LineWidth', 2, 'Color', [.8,0,0]);
title('Norm at each sample of Data vs. Minimum Norm solution.')