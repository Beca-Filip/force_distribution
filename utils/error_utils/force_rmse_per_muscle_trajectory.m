function [tensor_rmse] = force_rmse_per_muscle_trajectory(f1, f2, sample_list, trial_list, speed_list, leg_list)


% nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

tensor_rmse = zeros([size(f1, 1), ntrials, nspeeds, nlegs]);

% For all trials, speeds, and legs
for itr = 1 : ntrials
    for isp = 1 : nspeeds
        for ile = 1 : nlegs
            % For all muscles
            for mm = 1 : size(f1, 1)
                tensor_rmse(mm, itr, isp, ile) = rmse(f1(mm, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile)), ...
                                                      f2(mm, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile)));
            end
        end
    end
end


end