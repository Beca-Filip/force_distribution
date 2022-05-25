close all;
clear all;
clc;

% Bilevel results directory
bilevel_results_dir = '..\..\bilevel_optim_results\normalized_grid_search_5_partitions_0_60.mat';
load(bilevel_results_dir);

% Find the cost function parametrization producing least error
[~, min_ind] = min(E);
alpha_bilev = alpha(min_ind);

clearvars -except alpha_bilev

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);


% Exclude cost function 17
cf_exclude = 17;

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);
model.solver('ipopt');

% Perform IOC on these
trial_list = 1:10;
speed_list = 5;
leg_list = 1;
sample_list = 1:4:61;

% Active inequality thresh
aithresh = -1e-5;

% Create IO model
modelio = form_io_kktlsq_struct_normalized(data, model, vars, sample_list, trial_list, speed_list, leg_list);

% Process IO model
modelio = process_io_kkt_model_normalized(modelio, -1e-1);

% Analyze IO model
[~, ic_cf] = indexing_function(modelio, 'cf_parameters', [], [], [], []);
min(svd(modelio.regressor(:, ic_cf)))
min(svd(modelio.regressor))
[u,s,v] = svd(modelio.regressor);

solv.z = v(:, end);
[flag, solv.alpha] = positive_negative_orthant_check(normalize_vector(v(1:modelio.num_cf, end), 1));
solv.alpha.'

% Predict
solio = solve_io_kkt_model_normalized(modelio);

% DOC lists
doc_trial_list = 1:10;
doc_speed_list = 5;
doc_leg_list = 1;
doc_sample_list = 1:1:61;

% Redo-doc
Fout_KKT = DO_subroutine_normalized(solio.alpha, data, vars, model, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
% Redo-doc
Fout_bilevel = DO_subroutine_normalized(alpha_bilev, data, vars, model, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);

% Find means and standard deviations of data
[mean_f, std_f] = force_trajectory_statistics(data, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
stdpatch_f = cat(3, mean_f + 3 * std_f, flip(mean_f - 3 * std_f, 3));
stdpatch_t = cat(2, doc_sample_list-1, fliplr(doc_sample_list-1));

%% Find mean and standard deviations of KKT predictions
data_pred_KKT = data;
data_pred_KKT.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list) = Fout_KKT;
[mean_f_pred_KKT, std_f_pred_KKT] = force_trajectory_statistics(data_pred_KKT, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
stdpatch_f_pred_KKT = cat(3, mean_f_pred_KKT + 3 * std_f_pred_KKT, flip(mean_f_pred_KKT - 3 * std_f_pred_KKT, 3));
stdpatch_t_pred_KKT = cat(2, doc_sample_list-1, fliplr(doc_sample_list-1));

% Find mean and standard deviations of bilevel predictions
data_pred_bilev = data;
data_pred_bilev.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list) = Fout_bilevel;
[mean_f_pred_bilev, std_f_pred_bilev] = force_trajectory_statistics(data_pred_bilev, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
stdpatch_f_pred_bilev = cat(3, mean_f_pred_bilev + 3 * std_f_pred_bilev, flip(mean_f_pred_bilev - 3 * std_f_pred_bilev, 3));
stdpatch_t_pred_bilev = cat(2, doc_sample_list-1, fliplr(doc_sample_list-1));

%% Compare on plot
% Options for plotting
opts.title = @(n) {{sprintf('Muscle %d: RMSE$_{\\rm KKT}$ $= %.2f$N', n, ...
                    rmse(data.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred_KKT.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list)));...
                    sprintf('RMSE$_{\\rm Bilev}$ $= %.2f$N', ...
                    rmse(data.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred_bilev.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list)))},...
                    'interpreter', 'latex', 'FontSize', 15};
opts.xlabel = @(n) {'\% Gait Cycle', 'interpreter', 'latex', 'FontSize', 12};
opts.ylabel = @(n) {sprintf('$f_{%d}(t)$ [N]', n), 'interpreter', 'latex', 'FontSize', 13};
% Plot itself
figure('WindowState', 'Maximized');
% Data
patch_vector_quantities_opts(stdpatch_t, squeeze(stdpatch_f), [], [], [.2, .2, .5], 'FaceAlpha', .3);
plot_vector_quantities_opts(doc_sample_list-1, squeeze(mean_f), [], opts, 'b', 'LineWidth', 2);

patch_vector_quantities_opts(stdpatch_t_pred_KKT, squeeze(stdpatch_f_pred_KKT), [], [], [.5, .2, .2], 'FaceAlpha', .3);
plot_vector_quantities_opts(doc_sample_list-1, squeeze(mean_f_pred_KKT), [], [], 'r', 'LineWidth', 2);

patch_vector_quantities_opts(stdpatch_t_pred_bilev, squeeze(stdpatch_f_pred_bilev), [], [], [.2, .5, .2], 'FaceAlpha', .3);
plot_vector_quantities_opts(doc_sample_list-1, squeeze(mean_f_pred_bilev), [], [], 'g', 'LineWidth', 2);
% Overall title
sgtitle(sprintf('Overall RMSE$_{\\rm KKT}$ = %.4f [N], RMSE$_{\\rm Bilev}$ = %.4f [N]', ...
                rmse(data.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred_KKT.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list)), ...
                rmse(data.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred_bilev.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list)) ...
                ), 'interpreter', 'latex')
% Legend
lh = legend('$m_{\rm data} \pm 3 \sigma_{\rm data}$', '$m_{\rm data}$', '$m_{\rm KKT} \pm 3 \sigma_{\rm KKT}$', '$m_{\rm KKT}$', '$m_{\rm Bilev} \pm 3 \sigma_{\rm Bilev}$', '$m_{\rm Bilev}$', 'interpreter', 'latex', 'FontSize', 13);
title(lh, 'Legend');