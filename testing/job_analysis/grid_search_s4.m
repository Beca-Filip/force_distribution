close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Exclude cf
cf_exclude = [16, 17];

% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial_list = 1:10;
speed_list = 5;
leg_list = 2;
% sample_list = 1:4:61;
sample_list = 61:4:101;

% Get number of grid points
Ngrid = 3060;

% Calculate RMSE
tic
[E, alpha] = IO_grid_search_normalized(Ngrid, data, vars, model, sample_list, trial_list, speed_list, leg_list);
toc

prepresuffix = sprintf('Ngrid_%d', Ngrid);
presuffix = sprintf('leg_%d_speed_%d_', leg_list, speed_list);
suffix = datetimestr;

save(sprintf('..\\..\\opti_results\\job_patient4\\diffnorm_gs_%s_%s_%s.mat', prepresuffix, presuffix, suffix))