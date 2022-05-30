close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Create cost functions to exclude
% cf_exclude = [1, 5, 9, 16:17];
cf_exclude = [16:17];

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
alpha = zeros(17-length(cf_exclude), 1);
alpha(1:end) = 1;
alpha = normalize_vector(alpha, 1);

% Optimize with said vector
Fout = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% Set rng
rng(42);
% Fix noise levels
noise_lev = 1e-3;
noise_sigma = noise_lev * (data.fmax - data.fmin);
% Prealocate noisy Fout
Fout_noisy = zeros(size(Fout));
% For each muscle noise with different level
for ii = 1 : size(data.f, 1)
%     size(Fout(ii, :, :, :, :, :))
    Fout_noisy(ii, :, :, :, :, :) = Fout(ii, :, :, :, :, :) + noise_sigma(ii) * randn(size(Fout(ii, :, :, :, :, :)));
    mean(Fout_noisy(ii, :, :, :, :, :) - Fout(ii, :, :, :, :, :), 'all')
end
%%

% Create prediction structure
data_pred = data;
data_pred.f(:,:,sample_list,trial_list,speed_list,leg_list) = Fout_noisy;
% data_pred.f(:,:,sample_list,trial_list,speed_list,leg_list) = Fout;

%% Calculate RMSE
noisy_alpha = normalize_vector(alpha + 1e-1 * randn(size(alpha)), 1);
tic
[alpha_io, E] = IO_local_fmincon_search_normalized(noisy_alpha, data_pred, vars, model, sample_list, trial_list, speed_list, leg_list);
toc
