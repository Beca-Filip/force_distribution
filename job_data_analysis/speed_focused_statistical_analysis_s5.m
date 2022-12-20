close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\Optimization Model Data\Filtered_Patient5.mat';
load(data_dir);

% Perform analysis
% We use each of these at once
trial_list = 1:10;
sample_list = 1:101;
% We iterate over these
speed_list = [1 , 5];
leg_list = [1:2];

% Plot aesthetics
% Get a cmap according to what we want
cmap = linspecer(length(speed_list));
% linestyles
linstyl = {'--', ':', '-.', '-'};

% What to plot
% plot mean value flag
mean_flag = 1;
% plot standard deviation flag
std_flag = 1;

% Store results
% Size of the forces, but remove the dimension corresponding to trials
sizef = size(data.f);
sizef(4) = 1;
data_statistics.mean_f = zeros(sizef);
data_statistics.std_f = zeros(sizef);

% Iterate over speeds
for leg = leg_list
    
    % Plot options
%     opts.title = @(n) {sprintf('Muscle %d: $\\sigma = %.2f$N', n, mean(std_f(n, :, :))), 'interpreter', 'latex', 'FontSize', 15};
    opts.title = @(n) {sprintf('Muscle %d', n), 'interpreter', 'latex', 'FontSize', 15};
    opts.xlabel = @(n) {'\% Gait Cycle', 'interpreter', 'latex', 'FontSize', 12};
    opts.ylabel = @(n) {sprintf('$f_{%d}(t)$ [N]', n), 'interpreter', 'latex', 'FontSize', 13};
    % opts.ylim = @(n) {[min(data.f, [], 'all'), max(data.f, [], 'all')]};
    % opts.legend = @(n) {'$m \pm 3 \sigma$', '$m$', 'interpreter', 'latex'};
    
    % Create figure
    figure('Name', sprintf('Leg: %d', leg), 'Position', [0,0,1920,1080]);
    suptitle(sprintf("Leg: %d", leg));
    % Get for each leg calculate the statistics
    for speed = speed_list
        % Get statistics
        [mean_f, std_f] = force_trajectory_statistics(data, sample_list, trial_list, speed, leg);
        stdpatch_f = cat(3, mean_f + 3 * std_f, flip(mean_f - 3 * std_f, 3));
        stdpatch_t = cat(2, sample_list-1, fliplr(sample_list-1));

        % Store the results
        data_statistics.mean_f(:, :, :, 1, speed, leg) = mean_f;
        data_statistics.std_f(:, :, :, 1, speed, leg) = std_f;

        % Plot statistics
        
        % Indices for aesthetics
        cmap_index = find(speed == speed_list);
        linstyl_index = mod(find(speed == speed_list), length(linstyl));
        
        % Plot standard deviations only on request
        if std_flag
        % Darken the color for patch
        patch_vector_quantities_opts(stdpatch_t, squeeze(stdpatch_f), [], [], [.5, .5, .5], 'EdgeColor', cmap(cmap_index, :), 'FaceAlpha', .2, 'FaceColor', cmap(cmap_index, :), 'DisplayName', sprintf('$m_{speed_{%d}} \\pm 3\\sigma_{speed_{%d}}$', speed, speed));
        end
        % Plot the means only on request
        if mean_flag
        % Only get the opts once by plotting without opts until the last
        % pass
        if (speed ~= speed_list(end))
            [h_sp, ~] = plot_vector_quantities_opts(sample_list-1, squeeze(mean_f),  [], [], 'LineWidth', 2, ...
                'LineStyle', linstyl{linstyl_index}, 'Color', cmap(cmap_index, :), 'DisplayName', sprintf('$m_{speed_{%d}}$', speed));
        else
            [h_sp, ~] = plot_vector_quantities_opts(sample_list-1, squeeze(mean_f), [], opts, 'LineWidth', 2, ...
                'LineStyle', linstyl{linstyl_index}, 'Color', cmap(cmap_index, :), 'DisplayName', sprintf('$m_{speed_{%d}}$', speed));
        end
        end
        legen = legend('location', 'southeast');
        legen.Position(1) = h_sp{1, end}.Position(1);
        legen.Position(2) = h_sp{end, 1}.Position(2);
        legen.FontSize = 13;
        legen.Interpreter = 'Latex';
        
        if ~exist("graphs/", "dir")
           mkdir('graphs/') 
        end
        exportgraphics(gcf, sprintf("graphs/subject5_leg%d_speed%dvs%d.pdf", leg, speed_list),'contenttype','vector');
    end
end