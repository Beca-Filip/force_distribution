clear all;
close all;
clc;

% Data directory and loading
data_dir = '../../Optimization Model Data\Patient5.mat';
data = importdata(data_dir);

% Time vector
t = 0 : 100;
plotopts.xlabel = @(n) {ternary_operator(n >= 29, '\% Gait Cycle', ''), 'interpreter', 'latex'};
plotopts.ylabel = @(n) {sprintf("$f_{%d}$ [N]", n), 'interpreter', 'latex'};
% Plot all
for speed = 1 : size(data.f, 5)
% for speed = 1 : 1
    figure('WindowState', 'Maximized');
    % how many standard deviations to plot
    nstd = 1;

    mf1 = squeeze(mean(data.f(:, :, :, :, speed, 1), 4));
    sf1 = squeeze(std(data.f(:, :, :, :, speed, 1), 0, 4));

    mxf1 = squeeze(max(data.fmax(:, :, :, :, speed, 1), [], 4));
    mnf1 = squeeze(min(data.fmin(:, :, :, :, speed, 1), [], 4));

    plot_vector_quantities_opts_shape(t, mf1, [], plotopts, [], 'LineWidth', 2, 'Color', [1, 0, 0]);
    patch_vector_quantities_opts_shape([t, fliplr(t)], [mf1+nstd*sf1, fliplr(mf1-nstd*sf1)], [], [], [.7, .4, .4], [], 'FaceAlpha', .3, 'LineStyle', 'None');

    mf2 = squeeze(mean(data.f(:, :, :, :, speed, 2), 4));
    sf2 = squeeze(std(data.f(:, :, :, :, speed, 2), 0, 4));

    mxf2 = squeeze(max(data.fmax(:, :, :, :, speed, 2), [], 4));
    mnf2 = squeeze(min(data.fmin(:, :, :, :, speed, 2), [], 4));

    plot_vector_quantities_opts_shape(t, mf2, [], plotopts, [], 'LineWidth', 2, 'Color', [0, 0, 1]);
    patch_vector_quantities_opts_shape([t, fliplr(t)], [mf2+nstd*sf2, fliplr(mf2-nstd*sf2)], [], [], [.4, .4, .7], [], 'FaceAlpha', .3, 'LineStyle', 'None');

%         plot_vector_quantities_opts_shape(t, mxf1, [], [], [], '--', 'LineWidth', 2, 'Color', [.3, .3, 0]);
%         plot_vector_quantities_opts_shape(t, mnf1, [], [], [], '--', 'LineWidth', 2, 'Color', [0, .3, .3]);
    
    cmap1 = hot(size(data.f, 4) + 20);
    cmap1 = cmap1(11:end-10, :);
    cmap2 = cool(size(data.f, 4));
    for trial = 1 : size(data.f, 4)
        plot_vector_quantities_opts_shape(t, squeeze(data.f(:, :, :, trial, speed, 1)), [], [], [], ':', 'LineWidth', 0.5, 'Color', cmap1(trial, :));
        plot_vector_quantities_opts_shape(t, squeeze(data.f(:, :, :, trial, speed, 2)), [], [], [], ':', 'LineWidth', 0.5, 'Color', cmap2(trial, :));
    end
    
    for muscle = 1 : size(data.f, 1)
        % Which leg produces bigger muscle forces on average
        leg_force_diff = data.f(muscle, :, :, :, speed, 1) - data.f(muscle, :, :, :, speed, 2);
        avg_diff = mean(leg_force_diff, 'all');
        subplot(6, 6, muscle);
        box on;
        ax = gca;
        if avg_diff < 0
            ax.XColor = 'b';
            ax.YColor = 'b';
        else
            ax.XColor = 'r';
            ax.YColor = 'r';            
        end
    end
    
    figure('WindowState', 'Maximized');
    
    % Cross correlation parameters
    NumLags = 25;
    tlags = -NumLags : NumLags;
    % Prealocate structure for storing cross correlations
    leg_forces_xcorr = zeros(2*NumLags+1, size(data.f, 1), size(data.f, 4) * size(data.f, 4));
    for muscle = 1 : size(data.f, 1)
        for trial1 = 1 : size(data.f, 4)
            for trial2 = 1 : size(data.f, 4)
                % For every muscle and every trial calculate cross correlation
                [leg_forces_xcorr(:, muscle, (trial1 - 1) * size(data.f, 4) + trial2), ~] = crosscorr(squeeze(data.f(muscle, :, :, trial1, speed, 1)), squeeze(data.f(muscle, :, :, trial2, speed, 2)), 'NumLags', NumLags);
            end
        end
    end
    % Compute average cross-correlation across trials
    nstd = 1;
    avg_leg_forces_xcorr = mean(leg_forces_xcorr, 3);
    std_leg_forces_xcorr = std(leg_forces_xcorr, 0, 3);
    % Get maximum corr
    [max_leg_forces_xcorr, max_ind] = max(avg_leg_forces_xcorr);
    
    % Plot
    plotopts2.xlabel = @(n) {'Lag [\% Gait Cycle]', 'Interpreter', 'Latex'};
    plotopts2.ylabel = @(n) {sprintf("$C^{l_1, l_2}_{%d}$", n), 'interpreter', 'latex'};
    plotopts2.title = @(n) {sprintf("Max corr: $%.2f$; Lag: $%d$", max_leg_forces_xcorr(n), max_ind(n)), 'interpreter', 'latex'};
    plotopts2.yline = @(n) {0, 'k'};
    plotopts2.ylim = @(n) {[-1, 1]};
    plot_vector_quantities_opts_shape(tlags, avg_leg_forces_xcorr.', [], plotopts2, [], 'b');
    patch_vector_quantities_opts_shape([tlags, fliplr(tlags)], [avg_leg_forces_xcorr.' + std_leg_forces_xcorr.', fliplr(avg_leg_forces_xcorr.' - std_leg_forces_xcorr.')], [], [], [.5, .5, .5], [], 'FaceAlpha', .3, 'LineStyle', 'None');
    for muscle = 1 : size(data.f, 1)
        subplot(6, 6, muscle);
        if max_leg_forces_xcorr(muscle) > 0.85
            box on;
            ax = gca;
            ax.XColor = 'r';
            ax.YColor = 'r';
        elseif max_leg_forces_xcorr(muscle) > 0.65
            box on;
            ax = gca;
            ax.XColor = [0, .7, .3];
            ax.YColor = [0, .7, .3];
        end
    end
    
    
    figure('WindowState', 'Maximized');
    h = histogram(max_leg_forces_xcorr, linspace(0, 1, 21));
    xlabel('X-Corr. Coef.', 'Interpreter', 'latex');
    ylabel('\# Muscle Trajectories', 'Interpreter', 'latex');
    title('Histogram of average cross-correlation coefficients between trajectories of leg 1 and leg 2');
end