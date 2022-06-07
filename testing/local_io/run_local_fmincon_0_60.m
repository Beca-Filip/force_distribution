close all;
clear all;
clc;

% Data directory and loading
warning('off')
% name = 'diff_normalized_grid_search-5-partitions-0-60';
name = 'diff_normalized_grid_search-15-functions-5-partitions-0-60';
data_dir = strcat('..\..\bilevel_optim_results\diff_normalized_grid_search\', name, '.mat');
load(data_dir);
warning('on')


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


% Get sorted alphas
[Es, alphas] = IO_sort_rmses(E, alpha);

% Parameters of multiple local search 
n_run = 4;  % Points from which we'll run it
E_loc = zeros(n_run, 1);                    % Save results
alpha_loc = zeros(n_run, size(alpha, 2));   % Save results


% Run for the first 4
for ii = 1 : 4
    % Initialize
    E0 = Es(ii);
    alpha0 = alphas(ii, :);

    % Calculate RMSE
    tic
    [alpha_opt, E_opt] = IO_local_fmincon_search_normalized(alpha0, data, vars, model, sample_list, trial_list, speed_list, leg_list);
    toc
    
    % Store in array variables
    E_loc(ii) = E_opt;
    alpha_loc(ii, :) = alpha_opt;
end

% save_name = 'local-0-60';
save_name = 'local-15-functions-0-60';
save(strcat('..\..\bilevel_optim_results\diff_normalized_grid_search\', save_name, '.mat'));