function [tensor_ntcc] = force_temporal_correlation(f1, f2, sample_list, trial_list, speed_list, leg_list)


nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

tensor_ntcc = zeros([size(f1, [1, 2]), 1, ntrials, nspeeds, nlegs]);

% For all trials, speeds, and legs
for itr = 1 : ntrials
    for isp = 1 : nspeeds
        for ile = 1 : nlegs
            
            % Find the temporal mean of the signal
            mu1 = mean(f1(:, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile)), 3);
            mu2 = mean(f2(:, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile)), 3);
            % Find the temporal std of the signal
            sig1 = std(f1(:, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile)), 0, 3);
            sig2 = std(f2(:, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile)), 0, 3);
            
            % Find the temporal cross-correlation
            tcc = mean(...
                  (f1(:, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile)) - mu1) .* ...
                  (f2(:, :, sample_list, trial_list(itr), speed_list(isp), leg_list(ile)) - mu2),   ...
                  3);
              
            % Find the normalized temporal cross-correlation
            ntcc = (tcc ./ sig1) ./ sig2;
            
            % Store it
            tensor_ntcc(:,:,1,itr,isp,ile) = ntcc;
        end
    end
end


end