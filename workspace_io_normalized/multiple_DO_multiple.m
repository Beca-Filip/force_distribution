function [FOUT, SFLAG, EMSG] = multiple_DO_multiple(alpha_list, data, vars, model, sample_list, trial_list, speed_list, leg_list)

% Number of optimizations to do
nalphas = size(alpha_list, 1);
nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

% Create output structure
FOUT = zeros([size(vars.variables.f), nsamples, ntrials, nspeeds, nlegs, nalphas]);
% Prealocate success flag array
SFLAG = ones([size(vars.variables.f), nsamples, ntrials, nspeeds, nlegs, nalphas], 'logical');
% Prealocate error message string array
EMSG = strings([size(vars.variables.f), nsamples, ntrials, nspeeds, nlegs, nalphas]);

% Start alpha counter
for cntalphas = 1 : size(alpha_list, 1)
    
    % Assign alpha
    alpha = alpha_list(cntalphas, 1);
    
    % Do a mutliple DOs with same alpha
    [Fout, Sflag, Emsg] = DO_multiple(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

    % Store outputs
    FOUT(:, :, :, :, :, :, cntalphas) = Fout;
    SFLAG(:, :, :, :, :, :, cntalphas) = Sflag;
    EMSG(:, :, :, :, :, :, cntalphas) = Emsg;
end

end