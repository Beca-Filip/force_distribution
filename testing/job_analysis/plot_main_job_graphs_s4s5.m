close all;
clear all;
clc;

% Data directory and loading
data_dir_s1 = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir_s1);

% Modify normalization
data_s1 = data;
data_s1.J_min(:) = 0;
data_s1.J_max = data_s1.J_max ./ 1e3;

% Force bounds
data_s1.fmax(data_s1.fmax <= data_s1.f) = 1.003 * data_s1.f(data_s1.fmax <= data_s1.f);
data_s1.fmin(data_s1.fmin >= data_s1.f) = 0.997 * data_s1.f(data_s1.fmin >= data_s1.f);

epsil = 0.01;
data_s1.fmax(abs(data_s1.fmax - data_s1.f) < epsil) = data_s1.fmax(abs(data_s1.fmax - data_s1.fmin) < epsil) + epsil;
data_s1.fmin(abs(data_s1.f - data_s1.fmin) < epsil) = max(zeros(size(data_s1.fmin(abs(data_s1.f - data_s1.fmin) < epsil))), data_s1.fmin(abs(data_s1.f - data_s1.fmin) < epsil) - epsil);

% Exclude cf
cf_exclude = [16, 17];

% Create model
[model_s1, vars_s1] = form_casadi_model_normalized(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model_s1.solver('ipopt', sol_opt);

% Load IOC results
trials_s1 = load_job_data('ios4');

alphaArray_s1 = extract_properties_from_structs(trials_s1, 'alpha');
rmseArray_s1 = extract_properties_from_structs(trials_s1, 'err');

% Prealocate
predictedForces_s1 = nan(size(data_s1.f));
predictedRmses_s1 = nan([2, size(data_s1.f, 5:6)]);
predictedRmsesPerTrial_s1 = nan(size(data_s1.f, 4:6));

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
    alpha_curr = reshape(alphaArray_s1(leg_select, speed_select, phase_select, :, :), [], 1);

    % Compute solution
    fprintf("Computing solution phase: %d- speed: %d- leg: %d.\n", phase_select, speed_select, leg_select);
    Fout_curr = DO_subroutine_normalized(alpha_curr, data_s1, vars_s1, model_s1, sample_select, trial_select, speed_select, leg_select);

    % Store
    predictedForces_s1(:, :, sample_select, :, speed_select, leg_select) = Fout_curr;
    predictedRmses_s1(phase_select, speed_select, leg_select) = rmse(data_s1.f(:, :, sample_select, :, speed_select, leg_select), Fout_curr);
    
    % Check if rmse agree
    if abs(predictedRmses_s1(phase_select, speed_select, leg_select) - rmseArray_s1(leg_select, speed_select, phase_select)) > 1e-4
        fprintf("\tWARNING: Non-agreement of RMSEs.\n");
    end
end
end
end

for leg_select = [1, 2]
for speed_select = [1, 2, 3, 4, 5]
for individualTrialSelect = trial_select
    predictedRmsesPerTrial_s1(individualTrialSelect, speed_select, leg_select) = rmse(data_s1.f(:, :, :, individualTrialSelect, speed_select, leg_select), predictedForces_s1(:, :, :, individualTrialSelect, speed_select, leg_select));
end
end
end


% Data directory and loading
data_dir_s2 = '..\..\Optimization Model Data\Patient5.mat';
load(data_dir_s2);

% Modify normalization
data_s2 = data;
data_s2.J_min(:) = 0;
data_s2.J_max = data_s2.J_max ./ 1e3;

% Force bounds
data_s2.fmax(data_s2.fmax <= data_s2.f) = 1.003 * data_s2.f(data_s2.fmax <= data_s2.f);
data_s2.fmin(data_s2.fmin >= data_s2.f) = 0.997 * data_s2.f(data_s2.fmin >= data_s2.f);

epsil = 0.01;
data_s2.fmax(abs(data_s2.fmax - data_s2.f) < epsil) = data_s2.fmax(abs(data_s2.fmax - data_s2.fmin) < epsil) + epsil;
data_s2.fmin(abs(data_s2.f - data_s2.fmin) < epsil) = max(zeros(size(data_s2.fmin(abs(data_s2.f - data_s2.fmin) < epsil))), data_s2.fmin(abs(data_s2.f - data_s2.fmin) < epsil) - epsil);

% Exclude cf
cf_exclude = [16, 17];

% Create model
[model_s2, vars_s2] = form_casadi_model_normalized_s5(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model_s2.solver('ipopt', sol_opt);

% Load IOC results
trials_s2 = load_job_data('ios5');

alphaArray_s2 = extract_properties_from_structs(trials_s2, 'alpha');
rmseArray_s2 = extract_properties_from_structs(trials_s2, 'err');

% Prealocate
predictedForces_s2 = nan(size(data_s2.f));
predictedRmses_s2 = nan([2, size(data_s2.f, 5:6)]);
predictedRmsesPerTrial_s2 = nan(size(data_s2.f, 4:6));

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
    alpha_curr = reshape(alphaArray_s2(leg_select, speed_select, phase_select, :, :), [], 1);

    % Compute solution
    fprintf("Computing solution phase: %d- speed: %d- leg: %d.\n", phase_select, speed_select, leg_select);
    Fout_curr = DO_subroutine_normalized(alpha_curr, data_s2, vars_s2, model_s2, sample_select, trial_select, speed_select, leg_select);

    % Store
    predictedForces_s2(:, :, sample_select, :, speed_select, leg_select) = Fout_curr;
    predictedRmses_s2(phase_select, speed_select, leg_select) = rmse(data_s2.f(:, :, sample_select, :, speed_select, leg_select), Fout_curr);
    
    % Check if rmse agree
    if abs(predictedRmses_s2(phase_select, speed_select, leg_select) - rmseArray_s2(leg_select, speed_select, phase_select)) > 1e-4
        fprintf("\tWARNING: Non-agreement of RMSEs.\n");
    end
end
end
end


for leg_select = [1, 2]
for speed_select = [1, 2, 3, 4, 5]
for individualTrialSelect = trial_select
    predictedRmsesPerTrial_s2(individualTrialSelect, speed_select, leg_select) = rmse(data_s2.f(:, :, :, individualTrialSelect, speed_select, leg_select), predictedForces_s2(:, :, :, individualTrialSelect, speed_select, leg_select));
end
end
end


%% FIGURE: Weights

close all

% Objective functions
objective_fun_names = arrayfun(@(n) sprintf("$\\phi_{%d}$", n), 1:15, 'UniformOutput', false);
objective_fun_param_names = arrayfun(@(n) sprintf("$\\theta_{%d}$", n), 1:15, 'UniformOutput', false);
objective_fun_cmap = linspecer(15);

tick_locations = 1:15;
tick_strings = objective_fun_names;

% Speeds
speed_s1 = linspace(0.4, 0.8, 5);
speed_strings_s1 = arrayfun(@(n) sprintf("Speed=$$%.1f \\frac{\\textrm{m}}{\\textrm{s}}$$", speed_s1(n)), 1:5, 'UniformOutput', false);

speed_s2 = linspace(0.25, 0.65, 5);
speed_strings_s2 = arrayfun(@(n) sprintf("Speed=$$%.2f \\frac{\\textrm{m}}{\\textrm{s}}$$", speed_s2(n)), 1:5, 'UniformOutput', false);

% Phases
phase_strings = ["Stance phase.", "Swing phase."];

figure('WindowState', 'maximized');

% Subject 1, Leg 1, Phases 1 and 2
leg_select = 1;
speed_select = 2;
phase_select = 1;
ax1 = subplot(2, 4, 1);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s1(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
xlabel(speed_strings_s1(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\theta_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations);
xticklabels(tick_strings);
ax1.FontSize = 25; ax1.TickLabelInterpreter = 'latex';

leg_select = 1;
speed_select = 2;
phase_select = 2;
ax2 = subplot(2, 4, 2);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s1(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
xlabel(speed_strings_s1(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\theta_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations);
xticklabels(tick_strings);
ax2.FontSize = 25; ax2.TickLabelInterpreter = 'latex';


% Subject 1, Leg 2, Phases 1 and 2
leg_select = 2;
speed_select = 2;
phase_select = 1;
ax3 = subplot(2, 4, 3);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s1(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
xlabel(speed_strings_s1(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\theta_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations);
xticklabels(tick_strings);
ax3.FontSize = 25; ax3.TickLabelInterpreter = 'latex';

leg_select = 2;
speed_select = 2;
phase_select = 2;
ax4 = subplot(2, 4, 4);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s1(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
xlabel(speed_strings_s1(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\theta_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations);
xticklabels(tick_strings);
ax4.FontSize = 25; ax4.TickLabelInterpreter = 'latex';


% Subject 2, Leg 1, Phases 1 and 2
leg_select = 1;
speed_select = 3;
phase_select = 1;
ax5 = subplot(2, 4, 5);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s2(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
xlabel(speed_strings_s2(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\theta_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations);
xticklabels(tick_strings);
ax5.FontSize = 25; ax5.TickLabelInterpreter = 'latex';

leg_select = 1;
speed_select = 3;
phase_select = 2;
ax6 = subplot(2, 4, 6);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s2(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
xlabel(speed_strings_s2(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\theta_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations);
xticklabels(tick_strings);
ax6.FontSize = 25; ax6.TickLabelInterpreter = 'latex';


% Subject 2, Leg 2, Phases 1 and 2
leg_select = 2;
speed_select = 3;
phase_select = 1;
ax7 = subplot(2, 4, 7);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s2(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
xlabel(speed_strings_s2(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\theta_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations);
xticklabels(tick_strings);
ax7.FontSize = 25; ax7.TickLabelInterpreter = 'latex';

leg_select = 2;
speed_select = 3;
phase_select = 2;
ax8 = subplot(2, 4, 8);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s2(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
xlabel(speed_strings_s2(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\theta_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations);
xticklabels(tick_strings);
ax8.FontSize = 25; ax8.TickLabelInterpreter = 'latex';

% Legend
legend(objective_fun_param_names, 'Interpreter', 'latex', 'FontSize', 25, 'Units', 'normalized', 'Position', [0.9180 0.1100 0.0296 0.2839]);

% Annotations
tbSize = [0.1, 0.1];
xPosAnnotLeg1 = (ax1.Position(1) + ax2.Position(1) + ax2.Position(3) - tbSize(1)) / 2;
yPosAnnotLeg1 = (ax1.Position(2) + ax1.Position(4));
AnnotLeg1 = annotation('textbox', [xPosAnnotLeg1, yPosAnnotLeg1, tbSize], 'String', 'Non-paretic leg.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');

tbSize = [0.1, 0.1];
xPosAnnotLeg2 = (ax3.Position(1) + ax4.Position(1) + ax4.Position(3) - tbSize(1)) / 2;
yPosAnnotLeg2 = (ax3.Position(2) + ax3.Position(4));
AnnotLeg2 = annotation('textbox', [xPosAnnotLeg2, yPosAnnotLeg2, tbSize], 'String', 'Paretic leg.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');

tbSize = [0.1, 0.1];
xPosAnnotS1 = ax1.Position(1) - ax1.Position(3)/2 - tbSize(1)/2;
yPosAnnotS1 = (ax1.Position(2) + ax1.Position(4)/2 - tbSize(2)/2);
AnnotS1 = annotation('textbox', [xPosAnnotS1, yPosAnnotS1, tbSize], 'String', {'S1 : High'; 'functioning'; 'participant.'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');

tbSize = [0.1, 0.1];
xPosAnnotS2 = ax5.Position(1) - ax5.Position(3)/2 - tbSize(1)/2;
yPosAnnotS2 = (ax5.Position(2) + ax5.Position(4)/2 - tbSize(2)/2);
AnnotS2 = annotation('textbox', [xPosAnnotS2, yPosAnnotS2, tbSize], 'String', {'S2 : Low'; 'functioning'; 'participant.'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');

% % Separator lines
% vertPosBias = 0.025;
% xSeparatorS_1 = ax1.Position(1) - ax1.Position(3)/2;
% xSeparatorS_2 = ax4.Position(1) + 3*ax4.Position(3)/2;
% ySeparatorS = (ax1.Position(2) + ax5.Position(2) + ax5.Position(4)) / 2 - vertPosBias;
% annotation('line', [xSeparatorS_1 xSeparatorS_2], [ySeparatorS ySeparatorS], 'LineStyle', '--', 'LineWidth', 2, 'Color', [.5, .5, .5]);
% 
% horzPosBias = 0.01;
% xSeparatorL = (ax2.Position(1) + ax2.Position(3) + ax3.Position(1)) / 2 - horzPosBias;
% ySeparatorL_1 = ax2.Position(2) + ax2.Position(4);
% ySeparatorL_2 = ax6.Position(2);
% annotation('line', [xSeparatorL xSeparatorL], [ySeparatorL_1 ySeparatorL_2], 'LineStyle', '--', 'LineWidth', 2, 'Color', [.5, .5, .5]);

TileFigures


% Save figure
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/theta-vs-phase-leg-subj.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/theta-vs-phase-leg-subj.png'), 'ContentType', 'image', 'Resolution', 300);