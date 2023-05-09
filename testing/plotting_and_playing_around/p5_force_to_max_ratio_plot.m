clear all;
close all;
clc;

% Data directory and loading
data_dir = '../../Optimization Model Data\Patient5.mat';
data = importdata(data_dir);

% Time vector
t = 0 : 100;
muscle_indices = find(any(data.f > data.fmax, [2:6]));
plotopts.xlabel = @(n) {'\% Gait Cycle', 'interpreter', 'latex'};
plotopts.ylabel = @(n) {sprintf("$\\frac{f_{%d}}{{f_{%d}}_{\\rm max}}$", muscle_indices(n), muscle_indices(n)), 'interpreter', 'latex', 'fontsize', 20};
plotopts.title = @(n) {sprintf('$f_{%d}$', muscle_indices(n)), 'interpreter', 'latex'};
% Plot all
for speed = 1 : size(data.f, 5)
    for leg = 1 : size(data.f, 6)

        figure('WindowState', 'Maximized');
        hold on;
        
        sgtitle('Ratios of muscle force to maximum force: $\frac{f}{f_{\rm max}}$', 'interpreter', 'latex');
        
        % Plot limit of rations i.e. 1
        plot_vector_quantities_opts_shape(t([1, end]), ones(sum(any(data.f > data.fmax, [2:6])), 2), [], plotopts, [], 'LineWidth', 2, 'Color', [0, 0, 0]);
        
        % Modified maximum force estimates
        data.fmaxmod = max(data.f, data.fmax);
        
        for trial = 1 : size(data.f, 4)
            % Plot ratios
            plot_vector_quantities_opts_shape(t, squeeze(data.f(any(data.f > data.fmax, [2:6]), :, :, trial, speed, leg) ./ data.fmax(any(data.f > data.fmax, [2:6]), :, :, trial, speed, leg)), [], [], [], ':', 'LineWidth', 0.5);
            % Plot ratios with modified maximum estimates 
            plot_vector_quantities_opts_shape(t, squeeze(data.f(any(data.f > data.fmax, [2:6]), :, :, trial, speed, leg) ./ data.fmaxmod(any(data.f > data.fmax, [2:6]), :, :, trial, speed, leg)), [], [], [], '--', 'LineWidth', 0.5);
        end
    end
end