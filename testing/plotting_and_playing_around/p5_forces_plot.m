clear all;
close all;
clc;

% Data directory and loading
data_dir = '../../Optimization Model Data\Patient5.mat';
data = importdata(data_dir);

% Time vector
t = 0 : 100;
plotopts.xlabel = @(n) {ternary_operator(n >= 29, '\% Gait Cycle', ''), 'interpreter', 'latex'};
% Plot all
for speed = 1 : size(data.f, 5)
    for leg = 1 : size(data.f, 6)
        figure('WindowState', 'Maximized');
        mf = squeeze(mean(data.f(:, :, :, :, speed, leg), 4));
        sf = squeeze(std(data.f(:, :, :, :, speed, leg), 0, 4));
        
        mxf = squeeze(max(data.fmax(:, :, :, :, speed, leg), [], 4));
        mnf = squeeze(min(data.fmin(:, :, :, :, speed, leg), [], 4));
        
        plot_vector_quantities_opts_shape(t, mf, [], plotopts, [], 'LineWidth', 2, 'Color', [0, 0, 1]);
        patch_vector_quantities_opts_shape([t, fliplr(t)], [mf+3*sf, fliplr(mf-3*sf)], [], [], [.5, .5, .5], [], 'FaceAlpha', .3, 'LineStyle', 'None');
        
        plot_vector_quantities_opts_shape(t, mxf, [], [], [], '--', 'LineWidth', 2, 'Color', [.3, .3, 0]);
        plot_vector_quantities_opts_shape(t, mnf, [], [], [], '--', 'LineWidth', 2, 'Color', [0, .3, .3]);
        
        for trial = 1 : size(data.f, 4)
            plot_vector_quantities_opts_shape(t, squeeze(data.f(:, :, :, trial, speed, leg)), [], [], [], ':', 'LineWidth', 0.5);
        end
    end
end