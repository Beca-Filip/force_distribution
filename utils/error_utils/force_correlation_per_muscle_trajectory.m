function [tensor_pcc] = force_correlation_per_muscle_trajectory(f1, f2, sample_list, trial_list, speed_list, leg_list)


% nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

tensor_pcc = zeros([size(f1, 1), ntrials, nspeeds, nlegs]);

% For all trials, speeds, and legs
for itr = 1 : ntrials
    for isp = 1 : nspeeds
        for ile = 1 : nlegs
            % For all muscles
            for mm = 1 : size(f1, 1)
                tensor_pcc(mm, itr, isp, ile) = corr2(squeeze(f1(mm, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile))), ...
                                                      squeeze(f2(mm, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile))));
            end
        end
    end
end


end