function [mean_K, std_K] = stiffness_trajectory_statistics(data, sample_list, trial_list, speed_list, leg_list)
%STIFFNESS_TRAJECTORY_STATISTICS calculates the mean and standard deviation of
%force trajectories given specifications of which parts of the data to use.
%
%   [mean_f, std_f] = FORCE_TRAJECTORY_STATISTICS(data, sample_list, trial_list, speed_list, leg_list)

% Selected forces
sel_K = data.K(:, :, sample_list, trial_list, speed_list, leg_list);

% Create a mean along every dimension after the 3rd (meaning dont average
% across muscles nor samples)
mean_K = mean(sel_K, 4:length(size(sel_K)));

% Same for standard deviation
std_K = std(sel_K, 0, 4:length(size(sel_K)));
end