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
leg_list = 1;
sample_list = 1:4:61;

% Create weight vector
alpha = zeros(3, 1);
alpha(1:2) = 1;
alpha = normalize_vector(alpha, 1);

% Optimize with said vector
Fout = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);
% Create prediction structure
data_pred = data;
data_pred.f(:,:,sample_list,trial_list,speed_list,leg_list) = Fout;

%% Calculate RMSE
noisy_alpha = normalize_vector(alpha + 1e-2 * randn(size(alpha)), 1);
tic
[alpha_io, E] = IO_local_fmincon_search_normalized(noisy_alpha, data_pred, vars, model, sample_list, trial_list, speed_list, leg_list);
toc
