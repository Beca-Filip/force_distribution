%% FIGURE: Weights

close all

% Objective functions
% objective_fun_names = arrayfun(@(n) sprintf("$\\phi_{%d}$", n), 1:15, 'UniformOutput', false);
objective_fun_names = arrayfun(@(n) sprintf("$%d$", n), 1:2:15, 'UniformOutput', false);
objective_fun_param_names = arrayfun(@(n) sprintf("$\\omega_{%d}$", n), 1:15, 'UniformOutput', false);
objective_fun_cmap = linspecer(15);

tick_locations = 1:2:15;
tick_strings = objective_fun_names;

% Speeds
speed_s1 = linspace(0.4, 0.8, 5);
speed_strings_s1 = arrayfun(@(n) sprintf("Speed=$$%.1f \\frac{\\textrm{m}}{\\textrm{s}}$$", speed_s1(n)), 1:5, 'UniformOutput', false);

speed_s2 = linspace(0.25, 0.65, 5);
speed_strings_s2 = arrayfun(@(n) sprintf("Speed=$$%.2f \\frac{\\textrm{m}}{\\textrm{s}}$$", speed_s2(n)), 1:5, 'UniformOutput', false);

% Phases
phase_strings = ["Stance phase", "Swing phase"];

figure('WindowState', 'maximized');

% Subject 1, Leg 2, Phases 1 and 2
leg_select = 2;
speed_select = 2;
phase_select = 1;
ax1 = subplot(2, 4, 1);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s1(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = [0, 0, 0];
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
%xlabel(speed_strings_s1(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
xlabel("objective function", 'Interpreter', 'latex', 'FontSize', 25);
ylabel('inverse optimal control weight', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations); xtickangle(0);
xticklabels(tick_strings);
ax1.FontSize = 25; ax1.TickLabelInterpreter = 'latex';

leg_select = 2;
speed_select = 2;
phase_select = 2;
ax2 = subplot(2, 4, 2);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s1(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = [0, 0, 0];
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
%xlabel(speed_strings_s1(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
xlabel("objective function", 'Interpreter', 'latex', 'FontSize', 25);
ylabel('inverse optimal control weight', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations); xtickangle(0);
xticklabels(tick_strings);
ax2.FontSize = 25; ax2.TickLabelInterpreter = 'latex';

% Subject 2, Leg 2, Phases 1 and 2
leg_select = 2;
speed_select = 3;
phase_select = 1;
ax3 = subplot(2, 4, 5);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s2(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = [0, 0, 0];
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
%xlabel(speed_strings_s2(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
xlabel("objective function", 'Interpreter', 'latex', 'FontSize', 25);
ylabel('inverse optimal control weight', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations); xtickangle(0);
xticklabels(tick_strings);
ax3.FontSize = 25; ax3.TickLabelInterpreter = 'latex';

leg_select = 2;
speed_select = 3;
phase_select = 2;
ax4 = subplot(2, 4, 6);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, squeeze(alphaArray_s2(leg_select, speed_select, phase_select, ii, :)).');
    alphaBarPlot(ii).FaceColor = [0, 0, 0];
end
title(phase_strings(phase_select), 'FontSize', 25, 'Interpreter', 'latex');
%xlabel(speed_strings_s2(speed_select), 'Interpreter', 'latex', 'FontSize', 20);
xlabel("objective function", 'Interpreter', 'latex', 'FontSize', 25);
ylabel('inverse optimal control weight', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations); xtickangle(0);
xticklabels(tick_strings);
ax4.FontSize = 25; ax4.TickLabelInterpreter = 'latex';

TileFigures

% Annotations
tbSize = [0.1, 0.1];
% xPosAnnotS1 = ax1.Position(1) - ax1.Position(3)/2 - tbSize(1)/2;
xPosAnnotS1 = 0;
yPosAnnotS1 = (ax1.Position(2) + ax1.Position(4)/2 - tbSize(2)/2);
AnnotS1 = annotation('textbox', [xPosAnnotS1, yPosAnnotS1, tbSize], 'String', {'High'; 'functioning'; 'patient'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');

tbSize = [0.1, 0.1];
% xPosAnnotS2 = ax3.Position(1) - ax3.Position(3)/2 - tbSize(1)/2;
xPosAnnotS2 = 0;
yPosAnnotS2 = (ax3.Position(2) + ax3.Position(4)/2 - tbSize(2)/2);
AnnotS2 = annotation('textbox', [xPosAnnotS2, yPosAnnotS2, tbSize], 'String', {'Low'; 'functioning'; 'patient'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');


% Save figure
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/theta-vs-phase-leg-subj-isb-congress.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, sprintf('../../bilevel_optim_results/job_ioc_results/theta-vs-phase-leg-subj-isb-congress.png'), 'ContentType', 'image', 'Resolution', 300);