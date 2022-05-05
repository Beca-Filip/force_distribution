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
alpha(3) = 1;

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial = 9;
speed = 5;
leg = 1;
sample = 101;

% Determine active constraint
active = determine_active_inequalities_normalized(data, vars, model, sample, trial, speed, leg, 1e-6);

disp(active)