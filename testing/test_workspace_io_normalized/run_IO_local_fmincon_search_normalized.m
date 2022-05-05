close all;
clear all;
clc;

% Data directory and loading
warning('off')
name = 'normalized_grid_search-2022.05.04-09.24';
data_dir = strcat('..\..\opti_results\normalized_grid\', name, '.mat');
load(data_dir);
warning('on')

% Output directory and loading
out_dir = strcat('..\..\opti_results\normalized_local\', name);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Create model
[model, vars] = form_casadi_model_normalized();

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);

% Already defined in save
% % Perform IOC on these
% trial_list = 1:10;
% speed_list = 5;
% leg_list = 1;
% sample_list = 1:4:61;

%% Get sorted alphas
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
    OptDuration = toc

    % Name of current iteration and get path to save it
    name_curr = sprintf('iteration-%03d.mat', ii);
    save_path = strcat(out_dir, name_curr);
    
    % Save results
    save(save_path, 'E0', 'alpha0', 'E_opt', 'alpha_opt', 'OptDuration');
    
    % Store in array variables
    E_loc(ii) = E_opt;
    alpha_loc(ii, :) = alpha_opt;
end

% Save
name_curr = 'all-iterations.mat';
save_path = strcat(out_dir, name_curr);