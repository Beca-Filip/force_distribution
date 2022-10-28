close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\Optimization Model Data\Filtered_Patient4.mat';
load(data_dir);

% Perform analysis
% We use each of these at once
trial_list = 1:10;
sample_list = 1:101;
% We iterate over these
speed_list = [1 : 5];
leg_list = [1:2];

% Store results
% Size of the forces, but remove the dimension corresponding to trials
sizef = size(data.f);
sizef(4) = 1;
data_statistics.mean_f = zeros(sizef);
data_statistics.std_f = zeros(sizef);

% Iterate over speeds and legs
for leg = leg_list
for speed = speed_list
    % Get statistics
    [mean_f, std_f] = force_trajectory_statistics(data, sample_list, trial_list, speed, leg);
    stdpatch_f = cat(3, mean_f + 3 * std_f, flip(mean_f - 3 * std_f, 3));
    stdpatch_t = cat(2, sample_list-1, fliplr(sample_list-1));
    
    % Store the results
    data_statistics.mean_f(:, :, :, 1, speed, leg) = mean_f;
    data_statistics.std_f(:, :, :, 1, speed, leg) = std_f;
    
    % Plot options
    opts.title = @(n) {sprintf('Muscle %d: $\\sigma = %.2f$N', n, mean(std_f(n, :, :))), 'interpreter', 'latex', 'FontSize', 15};
    opts.xlabel = @(n) {'\% Gait Cycle', 'interpreter', 'latex', 'FontSize', 12};
    opts.ylabel = @(n) {sprintf('$f_{%d}(t)$ [N]', n), 'interpreter', 'latex', 'FontSize', 13};
    % opts.ylim = @(n) {[min(data.f, [], 'all'), max(data.f, [], 'all')]};
    % opts.legend = @(n) {'$m \pm 3 \sigma$', '$m$', 'interpreter', 'latex'};

    % Plot statistics
    % figure('WindowState', 'Maximized');
    figure('Name', sprintf('Leg: %d, Speed: %d', leg, speed), 'Position', [0,0,1920,1080]);
    patch_vector_quantities_opts(stdpatch_t, squeeze(stdpatch_f), [], [], [.5, .5, .5], 'FaceAlpha', .3, 'DisplayName', '$m \pm 3 \sigma$');
    plot_vector_quantities_opts(sample_list-1, squeeze(mean_f), [], opts, 'LineWidth', 2, 'DisplayName', '$m$');

    %
%     subplot;
%     leg = legend('interpreter', 'latex', 'FontSize', 13);
%     legend('interpreter', 'latex', 'FontSize', 13);
    legen = legend('location', 'southeast', '$m+3\sigma$','$m$');
    legen.Position(1) = legen.Position(1)+1/12;
    legen.FontSize = 13;
    legen.Interpreter = 'Latex';
%     title(leg, 'Legend');
end
end