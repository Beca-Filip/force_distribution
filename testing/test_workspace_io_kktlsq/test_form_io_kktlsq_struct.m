close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Exclude cost function 17
cf_exclude = 17;

% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);
model.solver('ipopt');

% Perform IOC on these
trial_list = 1:2;
speed_list = 5;
leg_list = 1;
sample_list = 1:2;

% Active inequality thresh
aithresh = -1e-5;

% Create IO model
modelio = form_io_kktlsq_struct_normalized(data, model, vars, sample_list, trial_list, speed_list, leg_list);