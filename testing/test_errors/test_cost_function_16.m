close all;
clear all;
clc;

rng(42);
% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Exclude cf
cf_exclude = [17];

% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 5;
% sol_opt.print_time =0;
% sol_opt.verbose = 1;
% sol_opt.ipopt.sb ='yes';
% sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
% sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);
% model.solver('ipopt');

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial_list = 1:10;
speed_list = 5;
leg_list = 1;
sample_list = 1:4:61;

% Set cost function weights
fdc = 4;
alpha = zeros(16, 1);
alpha(fdc) = 1;

% Perform DOC
Fout = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% Add noise to data which serves as initial conditions
sigma = 1e-1*(data.fmax(:,:,sample_list,trial_list,speed_list,leg_list) - data.fmin(:,:,sample_list,trial_list,speed_list,leg_list));
noise = sigma.*rand(size(data.f(:,:,sample_list,trial_list,speed_list,leg_list)));
data_noisy = data;
data_noisy.f(:,:,sample_list,trial_list,speed_list,leg_list) = data_noisy.f(:,:,sample_list,trial_list,speed_list,leg_list) + noise;

% Perform DOC with different (noisy) initial conditions
Fout_noisy = DO_subroutine_normalized(alpha, data_noisy, vars, model, sample_list, trial_list, speed_list, leg_list);

% Choose which trial to plot
randtrial = randi([1, 16], 1, 1);

%% Plot

% Max
max_ylim = max(max(Fout, Fout_noisy), [], [2 3 4 5 6]);
min_ylim = min(min(Fout, Fout_noisy), [], [2 3 4 5 6]);

opts.ylim = @(ii) {[min_ylim(ii), max_ylim(ii)]};
opts.xlabel = @(ii) {ternary_operator(ii >= 30, sprintf('\\%s Gait Cycle', '%'), ''), 'interpreter', 'latex'};
opts.ylabel = @(ii) {sprintf('$f_{%d}(t)$', ii), 'interpreter', 'latex'};

figure;
hold all;
plot_vector_quantities_opts(sample_list, squeeze(Fout(:, :, :, randtrial, :, :)), [], opts, 'LineWidth', 1);
plot_vector_quantities_opts(sample_list, squeeze(Fout_noisy(:, :, :, randtrial, :, :)), [], [], 'LineWidth', 1);
plot_vector_quantities_opts(sample_list, squeeze(data.fmin(:, :, sample_list, randtrial, speed_list, leg_list)), [], [], '--', 'LineWidth', 2, 'Color', [0, 0, 0]);
plot_vector_quantities_opts(sample_list, squeeze(data.fmax(:, :, sample_list, randtrial, speed_list, leg_list)), [], [], '--', 'LineWidth', 2, 'Color', [.5, .5, .5]);
legend(...
       sprintf('$f_0 = f_{\\rm human}$'), ...
       sprintf('$f_0 = f_{\\rm human} + {\\rm noise}$'), ...
       sprintf('$f_{\\rm min}$'), ...
       sprintf('$f_{\\rm max}$'), ...
       'interpreter', 'latex');
sgtitle(sprintf('Force trajectories obtained by minimizing $\\phi_{%d}$', fdc), 'interpreter', 'latex');
