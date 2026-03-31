function plotQuantilesAndPredictionsNoPhase(subjIdx, speedIdx, legIdx, data, predictedRmsePerTrial_1, predictedForces_1, predictedForces_2)
%PLOTQUANTILESANDPREDICTIONSNOPHASE  Plot force envelopes (min/max, IQR, median) and compare EMG vs IOC for the best trial.
%
% Inputs (expected shapes):
%   subjIdx, speedIdx, legIdx : indices
%   data.f            : forces, dims = [Nvar x ... x 101(time) x Ntrials x Nspeeds x Nlegs]
%   predictedRmsePerTrial(trial)
%   predictedForces    : same shape as data.f(:, :, :, :, sample, leg)
%
% Notes:
% - Uses custom helpers: patch_vector_quantities_opts_shape, plot_vector_quantities_opts_shape, ifelsef
% - X axis is gait cycle 0..100 (%)

% ---------- Select best (lowest RMSE) trial ----------
[~, bestTrialIdx] = min(squeeze(predictedRmsePerTrial_1));

% ---------- Aggregate envelopes & quantiles across trials ----------
% Take trial dimension (4th dim) min/max
forceMin = squeeze(min(data.f(:, :, :, :, speedIdx, legIdx), [], 4));
forceMax = squeeze(max(data.f(:, :, :, :, speedIdx, legIdx), [], 4));

% Explicit quartiles along trial dimension (Q1, median, Q3)
q = [0.25 0.50 0.75];
forceQuartiles = squeeze(quantile(data.f(:, :, :, :, speedIdx, legIdx), q, 4));
forceQ1 = squeeze(forceQuartiles(:, :, 1));
forceQ2 = squeeze(forceQuartiles(:, :, 2));  % median
forceQ3 = squeeze(forceQuartiles(:, :, 3));

% ---------- Axis/label options (per-subplot lambda signatures expected by your helpers) ----------
axisOpts = struct();
if subjIdx == 1
    % For participant 1, place xlabel on rows with n>=30 (per your subplot layout)
    axisOpts.xlabel = @(n) {ifelsef(n>=30, "\% Gait Cycle", ""), 'Interpreter','latex','FontSize',12};
else
    % For participant 2, place xlabel on rows with n>=29
    axisOpts.xlabel = @(n) {ifelsef(n>=29, "\% Gait Cycle", ""), 'Interpreter','latex','FontSize',12};
end
axisOpts.ylabel = @(n) {sprintf("$f^{%d}$ [N]", n), 'Interpreter','latex','FontSize',12};
axisOpts.xticks = @(n) {0:20:100};

% ---------- Plot ----------
x = 0:100;
xPatch = [x, fliplr(x)];

% Total range (min..max)
patch_vector_quantities_opts_shape( ...
    xPatch, [forceMin, fliplr(forceMax)], [], [], [.5 .5 .5], [6 6], ...
    'FaceAlpha', .2, 'LineStyle','none', ...
    'DisplayName', "$[f_{0\%}, f_{100\%}]$ [Total range]");

% Interquartile range (Q1..Q3)
patch_vector_quantities_opts_shape( ...
    xPatch, [forceQ1, fliplr(forceQ3)], [], [], [.2 .2 .8], [6 6], ...
    'FaceAlpha', .2, 'LineStyle','none', ...
    'DisplayName', "$[f_{25\%}, f_{75\%}]$ [IQR]");

% Median
plot_vector_quantities_opts_shape(x, forceQ2, [], [], [6 6], ...
    'LineWidth', 1, 'DisplayName', "$f_{50\%}$ [Median]");

% EMG (best trial)
plot_vector_quantities_opts_shape(x, squeeze(data.f(:, :, :, bestTrialIdx,  speedIdx, legIdx)), [], [], [6 6], ...
    'LineWidth', 2.5, 'LineStyle','--', "Color", [.2 .8 .1], ...
    'DisplayName', sprintf("$f_{\\rm EMG}$ [Gait cycle %d]", bestTrialIdx));

% IOC (best trial)
plot_vector_quantities_opts_shape(x, squeeze(predictedForces_1(:, :, :, bestTrialIdx)), [], axisOpts, [6 6], ...
    'LineWidth', 2.5, 'LineStyle',':', 'Color', [1 .435 1], ...
    'DisplayName', sprintf("$f_{\\rm IOC}$ [Gait cycle %d]", bestTrialIdx));

% Cross-Validation (best trial)
plot_vector_quantities_opts_shape(x, squeeze(predictedForces_2(:, :, :, bestTrialIdx)), [], axisOpts, [6 6], ...
    'LineWidth', 2.5, 'LineStyle','-.', 'Color', [0.86 0.08 0.24], ...
    'DisplayName', sprintf("$f_{\\rm CV}$ [Gait cycle %d]", bestTrialIdx));


% ---------- Legend ----------
legend('Position',[0.8315 0.0831 0.0742 0.1304], 'Interpreter','latex', 'FontSize',20);

% ---------- Axes cosmetics ----------
fig = gcf;
axList = findall(fig, 'Type','axes');
for k = 1:numel(axList)
    ax = axList(k);
    ax.FontSize = 15;
    ax.TickLabelInterpreter = 'latex';
end

% ---------- Title ----------
legLabels = ["Nonparetic","Paretic"];
sgtitle(sprintf("Participant S%d,  %s Leg", subjIdx, legLabels(legIdx)), ...
    'Interpreter','latex','FontSize',20);
end
