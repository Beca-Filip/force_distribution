function lin_ind = data_point_linear_index(modelio, sample, trial, speed, leg)
%DATA_POINT_LINEAR_INDEX determines the linear index of the data given its indices
%in the sample, trial, speed and leg space.
%
%   function lin_ind = DATA_POINT_LINEAR_INDEX(modelio, sample, trial, speed, leg)

    % Determine the linear index of the data point from its sample, trial
    % speed and leg index
    sample_ind = find(modelio.sample_list == sample);
    trial_ind = find(modelio.trial_list == trial);
    speed_ind = find(modelio.speed_list == speed);
    leg_ind = find(modelio.leg_list == leg);
    
    lin_ind = sub2ind([modelio.num_samples, modelio.num_trials, modelio.num_speeds, modelio.num_legs], ...
                       sample_ind, trial_ind, speed_ind, leg_ind);
end