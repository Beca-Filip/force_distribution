close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Create DOC model
[model, vars] = form_casadi_model_normalized();
model.solver('ipopt');

% Perform DOC-IOC on these
trial_list = 1:10;
speed_list = 1:5;
leg_list = 1;
sample_list = 1:101;

% Get statistics
[mean_f, std_f] = force_trajectory_statistics(data, sample_list, trial_list, speed_list, leg_list);
stdpatch_f = cat(3, mean_f + 3 * std_f, flip(mean_f - 3 * std_f, 3));
stdpatch_t = cat(2, sample_list-1, fliplr(sample_list-1));

% Plot options
opts.title = @(n) {sprintf('Muscle %d: $\\sigma = %.2f$', n, mean(std_f(n, :, :))), 'interpreter', 'latex'};
opts.xlabel = @(n) {'\% Gait Cycle', 'interpreter', 'latex'};
opts.ylabel = @(n) {sprintf('$f_{%d}(t)$', n), 'interpreter', 'latex', 'FontSize', 15};
% opts.ylim = @(n) {[min(data.f, [], 'all'), max(data.f, [], 'all')]};
opts.legend = @(n) {'$m \pm 3 \sigma$', '$m$', 'interpreter', 'latex'};

% Plot statistics
figure('WindowState', 'Maximized');
patch_vector_quantities_opts(stdpatch_t, squeeze(stdpatch_f), [], [], [.5, .5, .5], 'FaceAlpha', .3);
plot_vector_quantities_opts(sample_list-1, squeeze(mean_f), [], opts, 'LineWidth', 2);