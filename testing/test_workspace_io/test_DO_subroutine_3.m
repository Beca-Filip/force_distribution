close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Create model
[model, vars] = form_casadi_model();
model.solver('ipopt');

% Create weight vector
alpha = zeros(16, 1);
alpha(2) = 1;

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial_list = 1:10;
speed_list = 5;
leg_list = 1;
sample_list = 1:61;

% Calculate Forces
Fout = DO_subroutine(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

%% Plot
fig = figure;
fig.Position(1:2) = [100, 100];
fig.Position(3:4) = [1280 720];

% Which trial, speed, leg to plot
ttp = 3;
stp = 1;
ltp = 1;
% Time vector
Time = sample_list-1;
% Forces
plot_vector_quantities(Time, squeeze(data.f(:,:,sample_list,ttp,speed_list(stp),leg_list(ltp))), [], 'k', 'LineWidth', 2);
plot_vector_quantities(Time, squeeze(Fout(:,:,:,ttp)), [], 'r', 'LineWidth', 2);
% title
sgtitle(sprintf("RMSE: %.4f", rmse(squeeze(Fout(:,:,:,:)), squeeze(data.f(:,:,sample_list,:,speed_list(stp),leg_list(ltp)))) ));

exportgraphics(fig, '../../img/img_testing_workspace_io/test_DO_subroutine_3.png', 'ContentType', 'Image', 'Resolution', 400);