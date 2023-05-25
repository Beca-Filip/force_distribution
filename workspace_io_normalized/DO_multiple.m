function [Fout, Sflag, Emsg] = DO_multiple(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list)

% Number of optimizations to do
nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

% Create output structure
Fout = zeros([size(vars.variables.f), nsamples, ntrials, nspeeds, nlegs]);
% Prealocate success flag array
Sflag = ones([size(vars.variables.f), nsamples, ntrials, nspeeds, nlegs], 'logical');
% Prealocate error message string array
Emsg = strings([size(vars.variables.f, 1), size(vars.variables.f, 2), nsamples, ntrials, nspeeds, nlegs]);

% Set trial counter
cnttrials = 1;
% Loop over trials
for trial = trial_list
    
    % Reset speed counter
    cntspeeds = 1;
    % Loop over speeds
    for speed = speed_list
        
        % Reset leg counter
        cntlegs = 1;
        % Loop over legs
        for leg = leg_list
            
            % Reset sample counter
            cntsamples = 1;
            % Loop over samples
            for sample = sample_list
                
                % Do a single DO
                [fout, sflag, emsg] = DO_single(alpha, vars, model, data, sample, trial, speed, leg);
                
                % Store outputs
                Fout(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = fout;
                Sflag(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = sflag;
                Emsg(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = emsg;
                
                % Augment sample counter
                cntsamples = cntsamples + 1;
            end

            % Augment leg count
            cntlegs = cntlegs + 1;
        end
        
        % Augment speeds counter
        cntspeeds = cntspeeds + 1;
    end

    % Augment trial counter
    cnttrials = cnttrials + 1;
end

end