clear all;
close all;
clc;

% Raw data directory and loading
raw_data_dir = '..\Optimization Model Data\Filtered_Patient5.mat';
load(raw_data_dir);

% Exclude cfs
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

% Results dir
job_results_dir = "..\bilevel_optim_results\job_local_search\patient5";
results = job_load_results(job_results_dir);

data_ioc = make_predictions(data, results, model, vars);
results = compute_errors(data_ioc, results);
%%
generate_rmse_tables(results, "tables_s5/");
%% 
generate_alpha_tables(results, "tables_s5/");