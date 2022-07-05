
% Data copy
data_copy = data;

% Define lists
sample_list = 1:101;
trial_list = 1:10;
speed_list = 1:5;
leg_list = 1:2;

% Calculate stiffness
data_copy = compute_stiffness_wrap(data,sample_list,trial_list,speed_list,leg_list);

