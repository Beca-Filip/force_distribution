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
alpha = zeros(16, 1);
alpha(1) = 1;

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial_list = 5;
speed_list = 5;
leg_list = 1;
sample_list = 1:101;

% Calculate RMSE
E = IO_inner_loop_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);