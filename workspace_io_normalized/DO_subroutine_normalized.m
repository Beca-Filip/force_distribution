function Fout = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list)

% Number of optimizations to do
nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

% Create output structure
Fout = zeros([size(vars.variables.f), nsamples, ntrials, nspeeds, nlegs]);

% Counters
cntsamples = 1;
cnttrials = 1;
cntspeeds = 1;
cntlegs = 1;

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
            for k = sample_list

                % Set model parameters
                model = set_model_parameters(data, vars, model, k, trial, speed, leg);

                % Set model weights
                model = set_model_weights(alpha, vars, model);
                
                % Set cost function normalization
                model = set_model_normalization(data, vars, model);

                % Optimize
                sol = model.solve();

                % Take the solution forces
                f_opt = sol.value(vars.variables.f);

                % Store them in output
                Fout(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = f_opt;
                
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