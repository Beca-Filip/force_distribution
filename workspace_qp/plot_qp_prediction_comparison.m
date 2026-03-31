function plot_qp_prediction_comparison(subj_id, speed_idx, leg_idx, data, rmse_per_trial, f_predicted)
%PLOT_QP_PREDICTION_COMPARISON  Plot force envelopes and QP prediction for the best trial.
%
%   plot_qp_prediction_comparison(subj_id, speed_idx, leg_idx, data,
%                                 rmse_per_trial, f_predicted)
%
%   Draws one subplot per force component (6×6 grid) showing:
%     - Total range (min/max across trials) as a grey patch
%     - Interquartile range (Q1–Q3) as a blue patch
%     - Median across trials as a line
%     - Observed force for the best trial (lowest RMSE) as a dashed line
%     - QP-predicted force for the best trial as a dotted line
%
%   Inputs:
%     subj_id       - subject identifier (integer, used in title)
%     speed_idx     - speed index into data.f
%     leg_idx       - leg index into data.f
%     data          - struct with field:
%                       data.f  [n x 1 x nsamples x ntrials x nspeeds x nlegs]
%     rmse_per_trial - [ntrials x 1] per-trial RMSE of QP prediction
%     f_predicted   - [n x 1 x nsamples x ntrials] QP-predicted forces for
%                     this speed/leg (i.e. Fout(:,:,:,:,1,1) from QP_subroutine)

shape    = [6, 6];
nsamples = size(data.f, 3);
x        = linspace(0, 100, nsamples);
xpatch   = [x, fliplr(x)];

% ---- Select best (lowest RMSE) trial ------------------------------------
[~, best_trial_idx] = min(squeeze(rmse_per_trial));

% ---- Aggregate envelopes across trials ----------------------------------
f_sl = data.f(:, :, :, :, speed_idx, leg_idx);   % [n x 1 x nsamples x ntrials]

force_min = squeeze(min(f_sl, [], 4));             % [n x nsamples]
force_max = squeeze(max(f_sl, [], 4));             % [n x nsamples]

q = [0.25 0.50 0.75];
force_quartiles = squeeze(quantile(f_sl, q, 4));   % [n x nsamples x 3]
force_q1 = squeeze(force_quartiles(:, :, 1));      % [n x nsamples]
force_q2 = squeeze(force_quartiles(:, :, 2));      % [n x nsamples]  median
force_q3 = squeeze(force_quartiles(:, :, 3));      % [n x nsamples]

% ---- Axis options -------------------------------------------------------
last_row_start = (shape(1) - 1) * shape(2) + 1;   % first subplot index of bottom row

axis_opts.xlabel = @(curr) {ifelsef(curr >= last_row_start, "\% Gait Cycle", ""), ...
    'Interpreter', 'latex', 'FontSize', 12};
axis_opts.ylabel = @(curr) {sprintf("$f^{%d}$ [N]", curr), ...
    'Interpreter', 'latex', 'FontSize', 12};
axis_opts.xticks = @(~) {linspace(0, 100, 6)};

% ---- Plot ---------------------------------------------------------------

% Total range
patch_vector_quantities_opts_shape( ...
    xpatch, [force_min, fliplr(force_max)], [], [], [.5 .5 .5], shape, ...
    'FaceAlpha', .2, 'LineStyle', 'none', ...
    'DisplayName', "$[f_{0\%},\,f_{100\%}]$ [Total range]");

% IQR
patch_vector_quantities_opts_shape( ...
    xpatch, [force_q1, fliplr(force_q3)], [], [], [.2 .2 .8], shape, ...
    'FaceAlpha', .2, 'LineStyle', 'none', ...
    'DisplayName', "$[f_{25\%},\,f_{75\%}]$ [IQR]");

% Median
plot_vector_quantities_opts_shape(x, force_q2, [], [], shape, ...
    'LineWidth', 1, ...
    'DisplayName', "$f_{50\%}$ [Median]");

% Observed best trial
f_obs_best = squeeze(data.f(:, :, :, best_trial_idx, speed_idx, leg_idx));  % [n x nsamples]
plot_vector_quantities_opts_shape(x, f_obs_best, [], [], shape, ...
    'LineWidth', 2.5, 'LineStyle', '--', 'Color', [.2 .8 .1], ...
    'DisplayName', sprintf("$f_{\\rm obs}$ [Trial %d]", best_trial_idx));

% QP prediction best trial
f_pred_best = squeeze(f_predicted(:, :, :, best_trial_idx));                % [n x nsamples]
plot_vector_quantities_opts_shape(x, f_pred_best, [], axis_opts, shape, ...
    'LineWidth', 2.5, 'LineStyle', ':', 'Color', [1 .435 1], ...
    'DisplayName', sprintf("$f_{\\rm QP}$ [Trial %d]", best_trial_idx));

% ---- Legend and cosmetics -----------------------------------------------
legend('Position', [0.8315 0.0831 0.0742 0.1304], ...
    'Interpreter', 'latex', 'FontSize', 20);

fig = gcf;
for ax = findall(fig, 'Type', 'axes')'
    ax.FontSize = 15;
    ax.TickLabelInterpreter = 'latex';
end

leg_labels = ["Nonparetic", "Paretic"];
sgtitle(sprintf("Subject S%d,  Speed %d,  %s Leg", subj_id, speed_idx, leg_labels(leg_idx)), ...
    'Interpreter', 'latex', 'FontSize', 20);

end
