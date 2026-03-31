clear all;
close all;
clc;

% Data directory and loading
data_dir = '../../Optimization Model Data\Patient5.mat';
data = importdata(data_dir);

% Plot trial 1, speed 2, and leg 2 force before filtering
figure('WindowState', 'Maximized');
t = 0 : 100;
hold on;
plot_vector_quantities_opts_shape(t, squeeze(data.fmax(:,:,:,1,2,2)), [], [], []);

% Detect muscle indices that have the bug that F is bigger than FMAX
muscle_logical_indices = any(data.f > data.fmax, [2:6]);
muscle_indices = find(muscle_logical_indices);
% Filter the known error in patient 5 data for maximum forces
% Muscle 30, Sample 76/101, Trial 5, Speed 1, Leg 2
data.fmax(30, 1, 76, 5, 1, 2) = sum(data.fmax(30, 1, [75, 77], 5, 1, 2))/2;
% For all samples where FMAX is smaller than F you should replace FMAX with F
data.fmax(data.fmax < data.f) = data.fmax(data.fmax < data.f)  + ...
                                3 * (data.f(data.fmax < data.f) - data.fmax(data.fmax < data.f));

% Plot trial 1, speed 2, and leg 2 force before filtering
plot_vector_quantities_opts_shape(t, squeeze(data.fmax(:,:,:,1,2,2)), [], [], [], '--');
plot_vector_quantities_opts_shape(t, squeeze(data.f(:,:,:,1,2,2)), [], [], [], ':');

% Time vector
t = 0 : 100;
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
        plot_vector_quantities_opts_shape(t([1, end]), ones(sum(muscle_logical_indices), 2), [], plotopts, [], 'LineWidth', 2, 'Color', [0, 0, 0]);
                
        for trial = 1 : size(data.f, 4)
            % Plot ratios
            plot_vector_quantities_opts_shape(t, squeeze(data.f(muscle_logical_indices, :, :, trial, speed, leg) ./ data.fmax(muscle_logical_indices, :, :, trial, speed, leg)), [], [], [], ':', 'LineWidth', 0.5);
        end
    end
end

save('../../Optimization Model Data\Filtered_2_Patient5.mat', 'data');