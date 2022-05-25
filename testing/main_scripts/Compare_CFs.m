close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Exclude cost function 17
cf_exclude = 17;

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Create solver options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;

% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);
model.solver('ipopt', sol_opt);

% Perform DOC on these
doc_trial_list = 1:10;
doc_speed_list = 5;
doc_leg_list = 1;
% doc_sample_list = 1:4:61;
doc_sample_list = 61:4:101;

% Find means and standard deviations of data
[mean_f, std_f] = force_trajectory_statistics(data, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
stdpatch_f = cat(3, mean_f + 3 * std_f, flip(mean_f - 3 * std_f, 3));
stdpatch_t = cat(2, doc_sample_list-1, fliplr(doc_sample_list-1));

% Initialize data structure for predictions
data_pred = data;

% Predictions statistics
stat_pred = struct;
% Table of rmses
table_rmse_cc = table();

% Plot options
% Options for plotting
title_fun = @(n, f1, f2) {sprintf('Muscle %d: RMSE $= %.2f$N', n, rmse(f1, f2)), 'interpreter', 'latex', 'FontSize', 15};
%opts.title = @(n) title_fun(n,
%data.f(:,:,sample_list,trial_list,speed_list,leg_list),
%data_pred.f(:,:,sample_list,trial_list,speed_list,leg_list));
opts.xlabel = @(n) {'\% Gait Cycle', 'interpreter', 'latex', 'FontSize', 12};
opts.ylabel = @(n) {sprintf('$f_{%d}(t)$ [N]', n), 'interpreter', 'latex', 'FontSize', 13};

% For every cost function
n_cf = length(vars.functions.Jset);
fig = cell(n_cf, 1);
for cf_ii = 1 : n_cf
    
    % Put the corresponding weight to one
    alpha = zeros(n_cf, 1);
    alpha(cf_ii) = 1;
    
    % Do the DO routine
    Fout = DO_subroutine_normalized(alpha, data, vars, model, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
    
    % Stock the forces
    data_pred.f(:, :, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list) = Fout;
    
    % Get the statistics
    % Find mean and standard deviations of predictions
    [stat_pred(cf_ii).mean_f_pred, stat_pred(cf_ii).std_f_pred] = force_trajectory_statistics(data_pred, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
    stdpatch_f_pred = cat(3, stat_pred(cf_ii).mean_f_pred + 3 * stat_pred(cf_ii).std_f_pred, flip(stat_pred(cf_ii).mean_f_pred - 3 * stat_pred(cf_ii).std_f_pred, 3));
    stdpatch_t_pred = cat(2, doc_sample_list-1, fliplr(doc_sample_list-1));
    
    % Update title option
    opts.title = @(n) title_fun(n, data.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list));
    
    % Plot itself
    fig{cf_ii}=figure(cf_ii);
    fig{cf_ii}.WindowState = 'Maximized';
    % Data
    patch_vector_quantities_opts(stdpatch_t, squeeze(stdpatch_f), [], [], [.2, .2, .5], 'FaceAlpha', .3);
    plot_vector_quantities_opts(doc_sample_list-1, squeeze(mean_f), [], opts, 'b', 'LineWidth', 2);

    patch_vector_quantities_opts(stdpatch_t_pred, squeeze(stdpatch_f_pred), [], [], [.5, .2, .2], 'FaceAlpha', .3);
    plot_vector_quantities_opts(doc_sample_list-1, squeeze(stat_pred(cf_ii).mean_f_pred), [], [], 'r', 'LineWidth', 2);

    % Tabular RMSE
    table_rmse_cc.rmse(cf_ii) = rmse(data.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list));
    % Tabular CC
    table_rmse_cc.cc(cf_ii) = corr2(reshape(data.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), [], 1),...
                                      reshape(data_pred.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), [], 1));
    % Overall title
    sgtitle(sprintf('Overall RMSE = %.4f [N]', rmse(data.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list), data_pred.f(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list))), 'interpreter', 'latex')
    % Legend
    lh = legend('$m_{\rm data} \pm 3 \sigma_{\rm data}$', '$m_{\rm data}$', '$m_{\rm pred} \pm 3 \sigma_{\rm pred}$', '$m_{\rm pred}$', 'interpreter', 'latex', 'FontSize', 13);
    title(lh, 'Legend');
end

%% Save RMSE table - determine name from samples list
table_save_name = sprintf("rmse_cc_table_samples_%d-%d-%d", doc_sample_list(1)-1, mean(diff(doc_sample_list)), doc_sample_list(end)-1);
% Save RMSE table - fix table save directory
table_save_dir = sprintf("../../misc_results/misc_main_scripts/Compare_CFs");
% Get path
table_save_path = strcat(table_save_dir, "/", table_save_name, ".xlsx");

% Save table
writetable(table_rmse_cc, table_save_path);