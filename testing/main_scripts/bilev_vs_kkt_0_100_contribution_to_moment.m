clear all;
close all;
clc;


% Results of swing phase
sfr_dir = '../../bilevel_optim_results/bilev_kesh_comparison/local_bilev_vs_kesh_0_60.mat';
load(sfr_dir);

% Results of lift phase
lfr_dir = '../../bilevel_optim_results/bilev_kesh_comparison/local_bilev_vs_kesh_60_100.mat';
load(lfr_dir);

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Exclude cost function 17
cf_exclude = [16, 17];

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);
model.solver('ipopt');

% From the matrix of moment arms, calculate overall norm contribution
% factors
nA = multidim_column_norm(data.A, 2);
nAp = permute(nA, [2, 1, 3, 4, 5, 6]);

% Create KKT prediction
data_pred_KKT = data;
data_pred_KKT.f(:,:,swing_phase_results.indexing_list.samples, ...
                    swing_phase_results.indexing_list.trials, ...
                    swing_phase_results.indexing_list.speeds, ...
                    swing_phase_results.indexing_list.legs) = swing_phase_results.Fout_KKT;
data_pred_KKT.f(:,:,lift_phase_results.indexing_list.samples, ...
                    lift_phase_results.indexing_list.trials, ...
                    lift_phase_results.indexing_list.speeds, ...
                    lift_phase_results.indexing_list.legs) = lift_phase_results.Fout_KKT;

% Create bilevel prediction
data_pred_bilev = data;
data_pred_bilev.f(:,:,swing_phase_results.indexing_list.samples, ...
                      swing_phase_results.indexing_list.trials, ...
                      swing_phase_results.indexing_list.speeds, ...
                      swing_phase_results.indexing_list.legs) = swing_phase_results.Fout_bilevel;
data_pred_bilev.f(:,:,lift_phase_results.indexing_list.samples, ...
                      lift_phase_results.indexing_list.trials, ...
                      lift_phase_results.indexing_list.speeds, ...
                      lift_phase_results.indexing_list.legs) = lift_phase_results.Fout_bilevel;
                  
% Multiply all forces by the norm of the moment arms
data.f = data.f .* nAp;
data_pred_KKT.f = data_pred_KKT.f .* nAp;
data_pred_bilev.f = data_pred_bilev.f .* nAp;

% New indexing list
doc_trial_list = 1:10;
doc_speed_list = 5;
doc_leg_list = 1;
doc_sample_list = 1:4:101;

% Get statistics
[mean_f, std_f] = force_trajectory_statistics(data, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
stdpatch_f = cat(3, mean_f + 3 * std_f, flip(mean_f - 3 * std_f, 3));
stdpatch_t = cat(2, doc_sample_list-1, fliplr(doc_sample_list-1));

%% random trial for plot
for randtrial = 8

% Plot
% Options for plotting
% opts.title = @(n) {{sprintf('Muscle %d: RMSE$_{\\rm KKT}$ $= %.2f$N', n, ...
%                     rmse(data.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred_KKT.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list)));...
%                     sprintf('RMSE$_{\\rm Bilev}$ $= %.2f$N', ...
%                     rmse(data.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred_bilev.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list)))},...
%                     'interpreter', 'latex', 'FontSize', 12};
opts.xlabel = @(n) {ternary_operator(n>29, '\% Gait Cycle', []), 'interpreter', 'latex', 'FontSize', 12};
opts.ylabel = @(n) {sprintf('$\\| \\textbf{r}_{%d} \\|_2 \\cdot f_{%d}(t)$ [Nm]', n, n), 'interpreter', 'latex', 'FontSize', 12};
% opts.ylabel = @(n) {sprintf('$(\\| \\textbf{r} \\| \\cdot f)_{%d}(t)$ [Nm]', n), 'interpreter', 'latex', 'FontSize', 12};
opts.xline = @(n) {60, 'r--', 'LineWidth', 2, 'HandleVisibility', 'Off'};

% Plot itself
figure('WindowState', 'Maximized');
% Data
hold all;
plot_vector_quantities_opts(doc_sample_list, squeeze(data.f(:,:,doc_sample_list, randtrial, doc_speed_list, doc_leg_list)), [], opts, 'Color', [.0, .0, .0], 'LineWidth', 2, 'LineStyle', '--');
% plot_vector_quantities_opts(doc_sample_list, squeeze(data_pred_KKT.f(:,:,doc_sample_list, randtrial, doc_speed_list, doc_leg_list)), [], [], 'Color', [.2, .5, .2], 'LineWidth', 2);
plot_vector_quantities_opts(doc_sample_list, squeeze(data_pred_bilev.f(:,:,doc_sample_list, randtrial, doc_speed_list, doc_leg_list)), [], [], 'Color', [.2, .8, .2], 'LineWidth', 2 , 'LineStyle', '-.');
% Plot statistics
patch_vector_quantities_opts(stdpatch_t, squeeze(stdpatch_f), [], [], [.5, .5, .5], 'FaceAlpha', .3, 'LineStyle', 'None');
plot_vector_quantities_opts(doc_sample_list-1, squeeze(mean_f), [], [], 'LineWidth', 2, 'Color', [.2, .2, .8]);
% Plot separation of gait cycles

% Overall title
% sgtitle({...
%          sprintf('Overall RMSE$_{\\rm KKT}$ = %.4f [N], RMSE$_{\\rm Bilev}$ = %.4f [N]', ...
%                 rmse(data.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), data_pred_KKT.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list)), ...
%                 rmse(data.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), data_pred_bilev.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list)) ...
%                 ), ...
%          sprintf('Overall CC$_{\\rm KKT}$ = %.4f, CC$_{\\rm Bilev}$ = %.4f', ...
%                 corr2(reshape(data.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), [], 1), reshape(data_pred_KKT.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), [], 1)), ...
%                 corr2(reshape(data.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), [], 1), reshape(data_pred_bilev.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), [], 1)) ...
%                 )...
%         }, 'interpreter', 'latex')
% sgtitle({...
%          sprintf('Overall RMSE$_{\\rm Bilev}$ = %.4f [N] and', ...
%                  rmse(data.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), data_pred_bilev.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list)) ...
%                 ), ...
%          sprintf('CC$_{\\rm Bilev}$ = %.4f for gait cycle %d.', ...
%                 corr2(reshape(data.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), [], 1), reshape(data_pred_bilev.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), [], 1)), ...
%                 randtrial ...
%                 )...
%         }, 'interpreter', 'latex')
sgtitle({...
         sprintf('Overall RMSE$_{\\rm Bilev}$ = %.4f [Nm] and CC$_{\\rm Bilev}$ = %.4f for gait cycle %d.', ...
                 rmse(data.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), data_pred_bilev.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list)), ...
                corr2(reshape(data.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), [], 1), reshape(data_pred_bilev.f(:,:,doc_sample_list,randtrial,doc_speed_list,doc_leg_list), [], 1)), ...
                randtrial ...
                )...
        }, 'interpreter', 'latex')
% Legend
% lh = legend('$f_{\rm Human}(t)$', '$f_{\rm KKT}(t)$', '$f_{\rm Bilevel}(t)$', '$m \pm 3 \sigma$', '$m$', 'interpreter', 'latex', 'FontSize', 12);
lh = legend(sprintf('$\\| \\textbf{r}^{(%d)}\\| \\cdot f^{(%d)}_{\\rm Human}(t)$', randtrial, randtrial),...
            sprintf('$\\| \\textbf{r}^{(%d)}\\| \\cdot f^{(%d)}_{\\rm Bilevel}(t)$', randtrial, randtrial),...
            '$m_{\| \textbf{r}\| \cdot f_{\rm Human}} \pm 3 \sigma_{\| \textbf{r}\| \cdot f_{\rm Human}}$',...
            '$m_{\| \textbf{r}\| \cdot f_{\rm Human}}$',...
            'interpreter', 'latex', 'FontSize', 12);
title(lh, 'Legend');
end


%% Plot the alphas

% Alphas
bar1 = [swing_phase_results.alpha_bilev;swing_phase_results.solio.alpha.'].';
bar2 = [lift_phase_results.alpha_bilev;lift_phase_results.solio.alpha.'].';

% xlabels
xloc = 1 : length(swing_phase_results.alpha_bilev);
xname = arrayfun(@(n) sprintf('\\omega_{%d}', n), xloc, 'uniformoutput', false);

% Plot
figure('WindowState', 'Maximized');

% Swing phase
subplot(2, 1, 1);
hb1 = bar(bar1);
xticks(xloc);
xticklabels(xname);
legend('Bilevel', 'KKT', 'interpreter', 'latex', 'fontsize', 15);
title('Swing Phase', 'interpreter', 'latex', 'fontsize', 15);

% Lift phase
subplot(2, 1, 2);
hb2 = bar(bar2);
xticks(xloc);
xticklabels(xname);
legend('Bilevel', 'KKT', 'interpreter', 'latex', 'fontsize', 15);
title('Lift Phase', 'interpreter', 'latex', 'fontsize', 15);