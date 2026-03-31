function f_mean = compute_f_mean(f_data)
%COMPUTE_F_MEAN Compute the mean force vector across all conditions.
%
%   f_mean = COMPUTE_F_MEAN(f_data)
%
%   f_data is the force array with dimensions [n, 1, nsamples, ntrials,
%   nspeeds, nlegs].  All conditions are pooled so that f_mean represents
%   the global mean force across every sample, trial, speed and leg.
%
%   Output:
%     f_mean - (n x 1) mean force vector

n = size(f_data, 1);
F = reshape(f_data, n, []);     % [n x N], N = all conditions pooled
f_mean = mean(F, 2);            % [n x 1]

end
