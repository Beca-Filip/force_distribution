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
alpha(1) = 1;

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial_list = [1, 5];
speed_list = [2, 4];
leg_list = [1, 2];
sample_list = 1:2:101;

% Calculate Forces
Fout = DO_subroutine(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);