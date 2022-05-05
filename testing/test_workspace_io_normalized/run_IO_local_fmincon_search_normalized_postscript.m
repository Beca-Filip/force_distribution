sample_list_rerun = 1:61;
trial_list_rerun = 5;
speed_list_rerun = 5;
leg_list_rerun = 1;

% The output force
Fout = DO_subroutine_normalized(alpha_opt, data, vars, model, sample_list_rerun, trial_list_rerun, speed_list_rerun, leg_list_rerun)

% 
figure;
plot_vector_quantities(sample_list_rerun-1, squeeze(data.f(:,:,sample_list_rerun, trial_list_rerun, speed_list_rerun, leg_list_rerun)), ...
                [], 'LineWidth', 2, 'Color', 'k');
plot_vector_quantities(sample_list_rerun-1, squeeze(Fout), ...
                [], 'LineWidth', 2, 'Color', 'r');