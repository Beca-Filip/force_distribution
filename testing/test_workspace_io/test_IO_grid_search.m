close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Create model
[model, vars] = form_casadi_model();

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
trial_list = 5;
speed_list = 5;
leg_list = 1;
sample_list = 1:4:61;

% Get number of grid points
Ngrid = 17;

% Calculate RMSE
tic
[E, alpha] = IO_grid_search(Ngrid, data, vars, model, sample_list, trial_list, speed_list, leg_list);
toc