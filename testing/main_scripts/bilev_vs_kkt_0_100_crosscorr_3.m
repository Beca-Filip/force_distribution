% Find tensor of normalized temporal cross-correlation
tensor_ntcc = force_temporal_correlation(data.f, data_pred_bilev.f, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);

% Find mean across trials
m_trials_tensor_ntcc = mean(tensor_ntcc, 4);

% Plot
fig = figure;
fig.WindowState = 'Maximized';
bar3(m_trials_tensor_ntcc);

xlabel('Muscle axis');

ytickloc = 1:35;
yticklab = arrayfun(@(ii) sprintf('%d', ii), ytickloc, 'UniformOutput', false);
yticks(ytickloc);
yticklabels(yticklab);
yaxisproperties = get(gca, 'Xaxis');
yaxisproperties.TickLabelInterpreter = 'Latex';
ylim([0, 35])

% Find mean across all
m_all_tensor_ntcc = mean(tensor_ntcc, 'all');

fprintf('Mean of normalized temporal cross-correlation coefficient, averaged across all muscles and all trials is CC = %.4f.\n', m_all_tensor_ntcc);