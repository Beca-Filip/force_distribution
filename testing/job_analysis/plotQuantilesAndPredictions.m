function plotQuantilesAndPredictions(subjIdx, speedIdx, legIdx, data, predictedRmsePerTrial, predictedForces)

[~, minRmseTrialIdx] = min(squeeze(predictedRmsePerTrial(:, speedIdx, legIdx)));

lb_f = squeeze(min(data.f(:, :, :, :, speedIdx, legIdx), [], 4));
ub_f = squeeze(max(data.f(:, :, :, :, speedIdx, legIdx), [], 4));

quantiles_f = squeeze(quantile(data.f(:, :, :, :, speedIdx, legIdx), 3, 4));
one_quart_f = squeeze(quantiles_f(:, :, 1));
median_f = squeeze(quantiles_f(:, :, 2));
three_quart_f = squeeze(quantiles_f(:, :, 3));

% Participant 1 / Subject 4
if subjIdx == 1
    opts.xlabel = @(n) {ifelsef(n>=30, "\% Gait Cycle", ""), 'Interpreter', 'latex', 'FontSize', 12};
% Participant 2 / Subject 5
else
    opts.xlabel = @(n) {ifelsef(n>=29, "\% Gait Cycle", ""), 'Interpreter', 'latex', 'FontSize', 12};
end
opts.ylabel = @(n) {sprintf("$f^{%d}$ [N]", n), 'Interpreter', 'latex', 'FontSize', 12};
opts.xticks = @(n) {0:20:100};

% patch_vector_quantities_opts_shape([0:100, 100:-1:0], [lb_f, fliplr(ub_f)], [], [], [.5, .5, .5], [6, 6], 'FaceAlpha', .2, 'LineStyle', 'none', 'DisplayName', "$[f_{\\textrm{min}}, f_{\\textrm{max}}]$");
patch_vector_quantities_opts_shape([0:100, 100:-1:0], [lb_f, fliplr(ub_f)], [], [], [.5, .5, .5], [6, 6], 'FaceAlpha', .2, 'LineStyle', 'none', 'DisplayName', "$[f_{0\%}, f_{100\%}]$ [Total range region]");
patch_vector_quantities_opts_shape([0:100, 100:-1:0], [one_quart_f, fliplr(three_quart_f)], [], [], [.2, .2, .8], [6, 6], 'FaceAlpha', .2, 'LineStyle', 'none', 'DisplayName', "$[f_{25\%}, f_{75\%}]$ [Partial range region]");
plot_vector_quantities_opts_shape(0:100, median_f, [], [], [6, 6], 'LineWidth', 1, 'DisplayName', "$f_{50\%}$ [Median trajectory]");
plot_vector_quantities_opts_shape(0:100, squeeze(data.f(:, :, :, minRmseTrialIdx, speedIdx, legIdx)), [], [], [6, 6], 'LineWidth', 2.5, 'LineStyle', '--', 'DisplayName', sprintf("$f_{\\rm EMG}$ [Selected gait cycle %d]", minRmseTrialIdx));
plot_vector_quantities_opts_shape(0:100, squeeze(predictedForces(:, :, :, minRmseTrialIdx, speedIdx, legIdx)), [], opts, [6, 6], 'LineWidth', 2.5, 'LineStyle', ':', 'Color', [.2, .8, .1], 'DisplayName', sprintf("$f_{\\rm IOC}$ [Selected gait cycle %d]", minRmseTrialIdx));
% plot_vector_quantities_opts_shape(0:100, squeeze(data.f(:, :, :, minRmseTrialIdx, speedIdx, legIdx)), [], [], [6, 6], 'LineWidth', 2.5, 'LineStyle', '--', 'DisplayName', sprintf("$f_{\\rm EMG}^{(l=%d, s=%d, g=%d)}$", legIdx, speedIdx, minRmseTrialIdx));
% plot_vector_quantities_opts_shape(0:100, squeeze(predictedForces(:, :, :, minRmseTrialIdx, speedIdx, legIdx)), [], opts, [6, 6], 'LineWidth', 2.5, 'LineStyle', ':', 'Color', [.2, .8, .1], 'DisplayName', sprintf("$f_{\\rm IOC}^{(l=%d, s=%d, g=%d)}$", legIdx, speedIdx, minRmseTrialIdx));

% make legend
% legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 15);
legend('Position', [0.8315 0.0831 0.0742 0.1304], 'Interpreter', 'latex', 'FontSize', 20);

% Get the handle to the current figure
fig = gcf;

% Get all axes in the figure
axesHandles = findall(fig, 'Type', 'axes');

% Loop through each axes handle
for k = 1:length(axesHandles)
    % Get the current axes handle
    ax = axesHandles(k);
    
    ax.FontSize = 15;
    ax.TickLabelInterpreter = 'latex';
end

% Title
if subjIdx == 1
    speedsValues = linspace(0.4, 0.8, 5);
else
    speedsValues = linspace(0.25, 0.65, 5);
end

legLabel = ["Nonparetic", "Paretic"];
% sgtitle(sprintf("Participant S%d,  %s Leg,  Speed %.2f $\\frac{\\rm m}{\\rm s}$", subjIdx, legLabel(legIdx), speedsValues(speedIdx)), 'interpreter', 'latex', 'fontsize', 20);
sgtitle(sprintf("Participant S%d,  %s Leg", subjIdx, legLabel(legIdx)), 'interpreter', 'latex', 'fontsize', 20);


end