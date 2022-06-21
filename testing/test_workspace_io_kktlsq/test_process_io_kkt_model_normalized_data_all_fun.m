close all;
clear all;
clc;

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

% Perform IOC on these
trial_list = 1:10;
speed_list = 5;
leg_list = 1;
% sample_list = 1:4:61;
sample_list = 61:4:101;

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
Fout = DO_subroutine_normalized(solio.alpha, data, vars, model, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);

% Find means and standard deviations of data
[mean_f, std_f] = force_trajectory_statistics(data, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
stdpatch_f = cat(3, mean_f + 3 * std_f, flip(mean_f - 3 * std_f, 3));
stdpatch_t = cat(2, doc_sample_list-1, fliplr(doc_sample_list-1));

% Find mean and standard deviations of predictions
data_pred = data;
data_pred.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list) = Fout;
[mean_f_pred, std_f_pred] = force_trajectory_statistics(data_pred, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
stdpatch_f_pred = cat(3, mean_f_pred + 3 * std_f_pred, flip(mean_f_pred - 3 * std_f_pred, 3));
stdpatch_t_pred = cat(2, doc_sample_list-1, fliplr(doc_sample_list-1));

%% Compare on plot
% Options for plotting
opts.title = @(n) {sprintf('Muscle %d: RMSE $= %.2f$N', n, rmse(data.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list),...
                                                                   data_pred.f(n,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list))),...
                                                                   'interpreter', 'latex', 'FontSize', 15};
opts.xlabel = @(n) {'\% Gait Cycle', 'interpreter', 'latex', 'FontSize', 12};
opts.ylabel = @(n) {sprintf('$f_{%d}(t)$ [N]', n), 'interpreter', 'latex', 'FontSize', 13};
% Plot itself
figure('WindowState', 'Maximized');
% Data
patch_vector_quantities_opts(stdpatch_t, squeeze(stdpatch_f), [], [], [.2, .2, .5], 'FaceAlpha', .3);
plot_vector_quantities_opts(doc_sample_list-1, squeeze(mean_f), [], opts, 'b', 'LineWidth', 2);

patch_vector_quantities_opts(stdpatch_t_pred, squeeze(stdpatch_f_pred), [], [], [.5, .2, .2], 'FaceAlpha', .3);
plot_vector_quantities_opts(doc_sample_list-1, squeeze(mean_f_pred), [], [], 'r', 'LineWidth', 2);

% Overall title
sgtitle({...
         sprintf('Overall RMSE = %.4f [N]', rmse(data.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list)));...
         sprintf('Overall CC = %.4f', corr2(reshape(data.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), [], 1), reshape(data_pred.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), [], 1)))...
         }, 'interpreter', 'latex')
% Legend
lh = legend('$m_{\rm data} \pm 3 \sigma_{\rm data}$', '$m_{\rm data}$', '$m_{\rm pred} \pm 3 \sigma_{\rm pred}$', '$m_{\rm pred}$', 'interpreter', 'latex', 'FontSize', 13);
title(lh, 'Legend');