close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;



% Force bounds
data.fmax(data.fmax <= data.f) = 1.003 * data.f(data.fmax <= data.f);
data.fmin(data.fmin >= data.f) = 0.997 * data.f(data.fmin >= data.f);

epsil = 0.01;
data.fmax(abs(data.fmax - data.f) < epsil) = data.fmax(abs(data.fmax - data.fmin) < epsil) + epsil;
data.fmin(abs(data.f - data.fmin) < epsil) = max(zeros(size(data.fmin(abs(data.f - data.fmin) < epsil))), data.fmin(abs(data.f - data.fmin) < epsil) - epsil);

% Exclude cf
cf_exclude = [16, 17];

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

% Load IOC results
trials = load_job_data('ios4');

alphaArray = extract_properties_from_structs(trials, 'alpha');
rmseArray = extract_properties_from_structs(trials, 'err');

% Prealocate
predictedForces = nan(size(data.f));
predictedRmses = nan([2, size(data.f, 5:6)]);
predictedRmsesPerTrial = nan(size(data.f, 4:6));

% Do every leg, speed and phase
for leg_select = [1, 2]
for speed_select = [1, 2, 3, 4, 5]
for phase_select = [1, 2]
    if phase_select == 1
        sample_select = 1:61;
    else
        sample_select = 61:101;
    end
    trial_select = 1:10;
    
    % Extract the current alpha
    alpha_curr = reshape(alphaArray(leg_select, speed_select, phase_select, :, :), [], 1);

    % Compute solution
    fprintf("Computing solution phase: %d- speed: %d- leg: %d.\n", phase_select, speed_select, leg_select);
    Fout_curr = DO_subroutine_normalized(alpha_curr, data, vars, model, sample_select, trial_select, speed_select, leg_select);

    % Store
    predictedForces(:, :, sample_select, :, speed_select, leg_select) = Fout_curr;
    predictedRmses(phase_select, speed_select, leg_select) = rmse(data.f(:, :, sample_select, :, speed_select, leg_select), Fout_curr);
    
    % Check if rmse agree
    if abs(predictedRmses(phase_select, speed_select, leg_select) - rmseArray(leg_select, speed_select, phase_select)) > 1e-4
        fprintf("\tWARNING: Non-agreement of RMSEs.\n");
    end
end
end
end

for leg_select = [1, 2]
for speed_select = [1, 2, 3, 4, 5]
for individualTrialSelect = trial_select
    predictedRmsesPerTrial(individualTrialSelect, speed_select, leg_select) = rmse(data.f(:, :, :, individualTrialSelect, speed_select, leg_select), predictedForces(:, :, :, individualTrialSelect, speed_select, leg_select));
end
end
end

phaseStrings = ["Stance", "Swing"];
%% FIGURE: OBJECTIVE WEIGHTS
close all
objective_fun_names = arrayfun(@(n) sprintf("$\\phi_{%d}$", n), 1:15, 'UniformOutput', false);
objective_fun_cmap = linspecer(15);

separatorLocations = 1.5:1:4.5;

speeds = linspace(0.4, 0.8, 5);
speed_strings = arrayfun(@(n) sprintf("Speed=$$%.1f \\frac{\\textrm{m}}{\\textrm{s}}$$", speeds(n)), 1:5, 'UniformOutput', false);
tick_locations = 1:5;
tick_strings = speed_strings;

% Colors
% Define the RGB values for the shades of blue and red
lightBlue = [31, 119, 180] / 255; % Hex: #1f77b4
darkBlue = [0, 80, 158] / 255;    % Hex: #00509e
lightRed = [255, 64, 14] / 255;  % Hex: #ff7f0e
darkRed = [214, 39, 40] / 255;    % Hex: #d62728

% Stack the colors in an array
rmseLinePlotCmap = [
    lightBlue;
    darkBlue;
    lightRed;
    darkRed
];

legLineStyles = ["--", ":"];
cnt = 1;    % plot counter
for leg_select = [1, 2]
for phase_select = [1, 2]
    plot_string = sprintf("Leg %d : %s phase.", leg_select, phaseStrings(phase_select));

    figure('WindowState', 'maximized');
    hold all;
    % Left axis
    alphaBarPlot = bar(squeeze(alphaArray(leg_select, :, phase_select, :, :)));
    for ii = 1 : length(alphaBarPlot)
        alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
    end
    title(plot_string, 'FontSize', 15, 'Interpreter', 'latex');
    
    
    % Plot separators
    xline(separatorLocations, 'k--', 'HandleVisibility', 'off');
    
    ylabel('$\theta_i$ [mag.]', 'FontSize', 15, 'Interpreter', 'latex');
    ylim([0, 1]);
    xticks(tick_locations);
    xticklabels(tick_strings);
    ax = gca; ax.FontSize = 15; ax.TickLabelInterpreter = "latex";

    % Right axis
    yyaxis right
    rmseLinePlot = plot(squeeze(rmseArray(leg_select, :, phase_select)), 'LineStyle', legLineStyles(leg_select), 'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 25, 'Color', rmseLinePlotCmap(cnt, :), 'DisplayName', plot_string);
    ylabel('RMSE [N]', 'FontSize', 15, 'Interpreter', 'latex');
    ylim(expand_interval([min(rmseArray, [], 'all'), max(rmseArray, [], 'all')], 1.2));
    ax = gca; ax.FontSize = 15; ax.TickLabelInterpreter = "latex";

    % Relevant function selection
    relevantFuns = any(squeeze(alphaArray(leg_select, :, phase_select, :, :)) >= 1e-2, 1);
    % Legend only for relevant funs
    legend(alphaBarPlot(relevantFuns), objective_fun_names{relevantFuns}, 'Interpreter', 'latex', 'FontSize', 15, 'Units', 'normalized', 'Position', [0.0110 0.1072 0.0615 0.3260]);
    
    % Increment plot counter
    cnt = cnt + 1;
end
end

TileFigures


% Save figures
cnt = 1;    % plot counter
for leg_select = [1, 2]
for phase_select = [1, 2]
    exportgraphics(figure(cnt), sprintf('../../bilevel_optim_results/job_ioc_results/patient4/theta-vs-speed-leg-%d-phase-%d.pdf', leg_select, phase_select), 'ContentType', 'vector');
    exportgraphics(figure(cnt), sprintf('../../bilevel_optim_results/job_ioc_results/patient4/theta-vs-speed-leg-%d-phase-%d.png', leg_select, phase_select), 'ContentType', 'image', 'Resolution', 300);
    cnt = cnt + 1;
end
end

%% FIGURE: RMSE
close all

separatorLocations = 1.5:1:4.5;

speeds = linspace(0.4, 0.8, 5);
speed_strings = arrayfun(@(n) sprintf("Speed=$$%.1f \\frac{\\textrm{m}}{\\textrm{s}}$$", speeds(n)), 1:5, 'UniformOutput', false);
tick_locations = 1:5;
tick_strings = speed_strings;

% Colors
% Define the RGB values for the shades of blue and red
lightBlue = [31, 119, 180] / 255; % Hex: #1f77b4
darkBlue = [0, 80, 158] / 255;    % Hex: #00509e
lightRed = [255, 64, 14] / 255;  % Hex: #ff7f0e
darkRed = [214, 39, 40] / 255;    % Hex: #d62728

% Stack the colors in an array
rmseLinePlotCmap = [
    lightBlue;
    darkBlue;
    lightRed;
    darkRed
];

legLineStyles = ["--", ":"];
figure('WindowState', 'maximized');
hold all;
cnt = 1;    % plot counter
for leg_select = [1, 2]
for phase_select = [1, 2]

    plot_string = sprintf("Leg %d : %s phase.", leg_select, phaseStrings(phase_select));
    rmseLinePlot = plot(squeeze(rmseArray(leg_select, :, phase_select)), 'LineStyle', legLineStyles(leg_select), 'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 25, 'Color', rmseLinePlotCmap(cnt, :), 'DisplayName', plot_string);

    % Increment plot counter
    cnt = cnt + 1;
end
end

% Plot separators
xline(separatorLocations, 'k--', 'HandleVisibility', 'off');

% Aesthetics
ylabel('RMSE [N]', 'FontSize', 15, 'Interpreter', 'latex');
ylim(expand_interval([min(rmseArray, [], 'all'), max(rmseArray, [], 'all')], 1.2));
xticks(tick_locations);
xticklabels(tick_strings);
ax = gca; ax.FontSize = 15; ax.TickLabelInterpreter = "latex";
legend('Units', 'normalized', 'Position', [0.4468 0.1348 0.1455 0.0703], 'Interpreter', 'latex');

TileFigures

% Save figures
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/patient4/rmse-vs-speed.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/patient4/rmse-vs-speed.png'), 'ContentType', 'image', 'Resolution', 300);

%% FIGURE: Best Fit
close all;

figure('WindowState', 'maximized');
hold all;

plotQuantilesAndPredictions(1, 2, 1, data, predictedRmsesPerTrial, predictedForces);
TileFigures;

% Save figures
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/patient4/force-vs-time-leg-1.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/patient4/force-vs-time-leg-1.png'), 'ContentType', 'image', 'Resolution', 300);

%% 
close all;

figure('WindowState', 'maximized');
hold all;

plotQuantilesAndPredictions(1, 2, 2, data, predictedRmsesPerTrial, predictedForces);
TileFigures;

% Save figures
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/patient4/force-vs-time-leg-2.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/patient4/force-vs-time-leg-2.png'), 'ContentType', 'image', 'Resolution', 300);