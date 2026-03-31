close all

% Objective functions
objective_fun_names = arrayfun(@(n) sprintf("$%d$", n), 1:2:15, 'UniformOutput', false);
objective_fun_param_names = arrayfun(@(n) sprintf("$\\omega_{%d}$", n), 1:15, 'UniformOutput', false);
objective_fun_cmap = linspecer(15);

tick_locations = 1:2:15;
tick_strings = objective_fun_names;

figure('WindowState', 'maximized');

% Subject 1, Leg 1
iocData = importdata("./ios4/whole-cycle-speed-2-leg-1.mat");
alphaArray = iocData.alpha; fprintf("alpha: "); fprintf("\t%.2f", alphaArray); fprintf("\n")
ax1 = subplot(2, 2, 1);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, alphaArray(ii));
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
xlabel("$\phi_i$", 'Interpreter', 'latex', 'FontSize', 25);
ylabel('$\omega_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations); xtickangle(0);
xticklabels(tick_strings);
ax1.FontSize = 25; ax1.TickLabelInterpreter = 'latex';

% Subject 1, Leg 2
iocData = importdata("./ios4/whole-cycle-speed-2-leg-2.mat");
alphaArray = iocData.alpha; fprintf("alpha: "); fprintf("\t%.2f", alphaArray); fprintf("\n")
ax3 = subplot(2, 2, 2);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, alphaArray(ii));
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
xlabel("$\phi_i$", 'Interpreter', 'latex', 'FontSize', 25);
ylabel('$\omega_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations); xtickangle(0);
xticklabels(tick_strings);
ax3.FontSize = 25; ax3.TickLabelInterpreter = 'latex';

% Subject 2, Leg 1
iocData = importdata("./ios5/whole-cycle-speed-3-leg-1.mat");
alphaArray = iocData.alpha; fprintf("alpha: "); fprintf("\t%.2f", alphaArray); fprintf("\n")
ax5 = subplot(2, 2, 3);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, alphaArray(ii));
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
xlabel("$\phi_i$", 'Interpreter', 'latex', 'FontSize', 25);
ylabel('$\omega_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations); xtickangle(0);
xticklabels(tick_strings);
ax5.FontSize = 25; ax5.TickLabelInterpreter = 'latex';

% Subject 2, Leg 2
iocData = importdata("./ios5/whole-cycle-speed-3-leg-2.mat");
alphaArray = iocData.alpha; fprintf("alpha: "); fprintf("\t%.2f", alphaArray); fprintf("\n")
ax7 = subplot(2, 2, 4);
hold all;
alphaBarPlot = gobjects(1, 15);
for ii = 1 : 15
    alphaBarPlot(ii) = bar(ii, alphaArray(ii));
    alphaBarPlot(ii).FaceColor = objective_fun_cmap(ii, :);
end
xlabel("$\phi_i$", 'Interpreter', 'latex', 'FontSize', 25);
ylabel('$\omega_i$ [mag.]', 'FontSize', 25, 'Interpreter', 'latex');
ylim([0, 1]);
xticks(tick_locations); xtickangle(0);
xticklabels(tick_strings);
ax7.FontSize = 25; ax7.TickLabelInterpreter = 'latex';

% Legend
legend(objective_fun_param_names, 'Interpreter', 'latex', 'FontSize', 25, 'Units', 'normalized', 'Position', [0.9180 0.1200 0.0296 0.2839]);

% % Annotations
% tbSize = [0.1, 0.1];
% xPosAnnotLeg1 = (ax1.Position(1) + ax2.Position(1) + ax2.Position(3) - tbSize(1)) / 2;
% yPosAnnotLeg1 = (ax1.Position(2) + ax1.Position(4));
% AnnotLeg1 = annotation('textbox', [xPosAnnotLeg1, yPosAnnotLeg1, tbSize], 'String', 'Non-paretic leg', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');
% 
% tbSize = [0.1, 0.1];
% xPosAnnotLeg2 = (ax3.Position(1) + ax4.Position(1) + ax4.Position(3) - tbSize(1)) / 2;
% yPosAnnotLeg2 = (ax3.Position(2) + ax3.Position(4));
% AnnotLeg2 = annotation('textbox', [xPosAnnotLeg2, yPosAnnotLeg2, tbSize], 'String', 'Paretic leg', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');
% 
% tbSize = [0.1, 0.1];
% xPosAnnotS1 = ax1.Position(1) - ax1.Position(3)/2 - tbSize(1)/2;
% yPosAnnotS1 = (ax1.Position(2) + ax1.Position(4)/2 - tbSize(2)/2);
% AnnotS1 = annotation('textbox', [xPosAnnotS1, yPosAnnotS1, tbSize], 'String', {'High'; 'functioning'; 'participant'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');
% 
% tbSize = [0.1, 0.1];
% xPosAnnotS2 = ax5.Position(1) - ax5.Position(3)/2 - tbSize(1)/2;
% yPosAnnotS2 = (ax5.Position(2) + ax5.Position(4)/2 - tbSize(2)/2);
% AnnotS2 = annotation('textbox', [xPosAnnotS2, yPosAnnotS2, tbSize], 'String', {'Low'; 'functioning'; 'participant'}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none', 'FontSize', 30, 'Interpreter', 'latex');
% 
% % % Separator lines
% % vertPosBias = 0.025;
% % xSeparatorS_1 = ax1.Position(1) - ax1.Position(3)/2;
% % xSeparatorS_2 = ax4.Position(1) + 3*ax4.Position(3)/2;
% % ySeparatorS = (ax1.Position(2) + ax5.Position(2) + ax5.Position(4)) / 2 - vertPosBias;
% % annotation('line', [xSeparatorS_1 xSeparatorS_2], [ySeparatorS ySeparatorS], 'LineStyle', '--', 'LineWidth', 2, 'Color', [.5, .5, .5]);
% % 
% % horzPosBias = 0.01;
% % xSeparatorL = (ax2.Position(1) + ax2.Position(3) + ax3.Position(1)) / 2 - horzPosBias;
% % ySeparatorL_1 = ax2.Position(2) + ax2.Position(4);
% % ySeparatorL_2 = ax6.Position(2);
% % annotation('line', [xSeparatorL xSeparatorL], [ySeparatorL_1 ySeparatorL_2], 'LineStyle', '--', 'LineWidth', 2, 'Color', [.5, .5, .5]);
% 
% TileFigures
% 
% % Rectangles
% annotation(gcf, 'rectangle', [0.01, 0.5125, .95, .4875], 'LineWidth', 2);
% annotation(gcf, 'rectangle', [0.01, 0.025, .95, .4875], 'LineWidth', 2);
% annotation(gcf, 'rectangle', [0.01, 0.5125, .08125, .4875], 'LineWidth', 2);
% annotation(gcf, 'rectangle', [0.01, 0.025, .08125, .4875], 'LineWidth', 2);
% annotation(gcf, 'rectangle', [0.01, 0.5125, .49, .4875], 'LineWidth', 2);
% annotation(gcf, 'rectangle', [0.01, 0.025, .49, .4875], 'LineWidth', 2);

% --- Final layout + annotations for 2x2 (one graph per participant & leg) ---

% Convenience handles
axTL = ax1;   % Subject 1, Leg 1 (top-left)
axTR = ax3;   % Subject 1, Leg 2 (top-right)
axBL = ax5;   % Subject 2, Leg 1 (bottom-left)
axBR = ax7;   % Subject 2, Leg 2 (bottom-right)
axs  = [axTL, axTR, axBL, axBR];

% Figure-wide bounds of axes block
% xLeft   = min(arrayfun(@(h) h.Position(1),               axs));
% xRight  = max(arrayfun(@(h) h.Position(1)+h.Position(3), axs));
% yBottom = min(arrayfun(@(h) h.Position(2),               axs));
% yTop    = max(arrayfun(@(h) h.Position(2)+h.Position(4), axs));
xLeft   = 0.1;
xRight  = 0.9;
yBottom = 0.05;
yTop    = 0.95;

% Column geometry (from top row)
col1_xc   = axTL.Position(1) + axTL.Position(3)/2;
col2_xc   = axTR.Position(1) + axTR.Position(3)/2;
col_gap_x = axTR.Position(1) - (axTL.Position(1)+axTL.Position(3));  %#ok<NASGU> % (not used, but handy)

% Row geometry
row1_y_top    = axTL.Position(2) + axTL.Position(4);
row1_y_bottom = axTL.Position(2);
row2_y_top    = axBL.Position(2) + axBL.Position(4);
row2_y_bottom = axBL.Position(2);

% ===== Column headers (legs) =====
tbSizeCol = [0.18, 0.06];  % [w, h]
yColTitle = row1_y_top + 0.01;  % small margin above top row
annotation('textbox', [col1_xc - tbSizeCol(1)/2, yColTitle, tbSizeCol], ...
    'String', 'Non-paretic leg', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'EdgeColor','none', ...
    'FontSize', 30, 'Interpreter','latex');

annotation('textbox', [col2_xc - tbSizeCol(1)/2, yColTitle, tbSizeCol], ...
    'String', 'Paretic leg', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'EdgeColor','none', ...
    'FontSize', 30, 'Interpreter','latex');

% ===== Row side labels (participants) =====
tbSizeRow = [0.12, 0.12];
xSide = axTL.Position(1) - tbSizeRow(1)/2;  % left of left column

annotation('textbox', [xSide, row1_y_bottom + (row1_y_top-row1_y_bottom)/2 - tbSizeRow(2)/2, tbSizeRow], ...
    'String', {'High','functioning','participant'}, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'EdgeColor','none', 'FontSize', 30, 'Interpreter','latex', 'Rotation', 90);

annotation('textbox', [xSide, row2_y_bottom + (row2_y_top-row2_y_bottom)/2 - tbSizeRow(2)/2, tbSizeRow], ...
    'String', {'Low','functioning','participant'}, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'EdgeColor','none', 'FontSize', 30, 'Interpreter','latex', 'Rotation', 90);

% ===== Separator lines =====
% Horizontal between rows
ySep = (row1_y_bottom + row2_y_top) / 2 - 0.025;
annotation('line', [xLeft, xRight], [ySep, ySep], ...
    'LineStyle','--', 'LineWidth',2, 'Color',[.5 .5 .5]);

% Vertical between columns
xSep = (axTL.Position(1)+axTL.Position(3) + axTR.Position(1)) / 2;
annotation('line', [xSep, xSep], [yBottom, yTop], ...
    'LineStyle','--', 'LineWidth',2, 'Color',[.5 .5 .5]);

TileFigures

exportgraphics(gcf, "..\..\bilevel_optim_results\job2_no_phase\s4s5-retrieved-cf.pdf", "ContentType", "vector");
exportgraphics(gcf, "..\..\bilevel_optim_results\job2_no_phase\s4s5-retrieved-cf.png", "ContentType", "image", "Resolution", 300);