close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Create cost functions to exclude
cf_exclude = [1, 3:5, 7:9, 11:17];

% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);
model.solver('ipopt');

% Create weight vector
alpha = zeros(3, 1);
alpha(1) = 1;
alpha = normalize_vector(alpha, 1);

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial_list = 1:10;
speed_list = 5;
leg_list = 1;
sample_list = 1:4:61;

% Get all the outputs
Fout = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

%% Calculate RMSE

% Store in structure
data_pred = data;
data_pred.f(:,:,sample_list,trial_list, speed_list, leg_list) = Fout;

% Do one iteration of IOC with right weights
E0 = IO_inner_loop_normalized(alpha, data_pred, vars, model, sample_list, trial_list, speed_list, leg_list);

%% Calculate RMSE

% Do one iteration of IOC with wrong weights
E1 = IO_inner_loop_normalized(normalize_vector(alpha + rand(size(alpha)), 1), data_pred, vars, model, sample_list, trial_list, speed_list, leg_list);