close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Filtered_2_Patient5.mat';
load(data_dir);

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Exclude cf
cf_exclude = [16, 17];

% Create model
[model, vars] = form_casadi_model_normalized_s5(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);

%% Grid

% Define parameter lists for which to optimize
alpha_list = prob_simplex_ndim(14, 7);
sample_list = 1 : 4 : 101;
trial_list = 1 : 10;
speed_list = 1 : 2;
leg_list = 1 : 2;

% Test
[FOUT, SFLAG, EMSG] = multiple_DO_multiple(alpha_list, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% Save
save('P5.mat', 'FOUT', 'SFLAG', 'EMSG');