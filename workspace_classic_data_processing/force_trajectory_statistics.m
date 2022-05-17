function [mean_f, std_f] = force_trajectory_statistics(data, sample_list, trial_list, speed_list, leg_list)
%FORCE_TRAJECTORY_STATISTICS calculates the mean and standard deviation of
%force trajectories given specifications of which parts of the data to use.
%
%   [mean_f, std_f] = FORCE_TRAJECTORY_STATISTICS(data, sample_list, trial_list, speed_list, leg_list)

% Selected forces
sel_f = data.f(:, :, sample_list, trial_list, speed_list, leg_list);

% Create a mean along every dimension after the 3rd (meaning dont average
% across muscles nor samples)
mean_f = mean(sel_f, 4:length(size(sel_f)));

% Same for standard deviation
std_f = std(sel_f, 0, 4:length(size(sel_f)));
end