function Fout = DO_subroutine_normalized_variable_phase(alpha, data, vars, model, sample_cell_list, trial_list, speed_list, leg_list)

% Number of trajectories to compute
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

% Each trial has a different number of samples so store outputs in cell
nsamples = zeros(ntrials, nspeeds, nlegs);
Fout = cell(ntrials, nspeeds, nlegs);
% For every trial
for trial = 1 : ntrials
    for speed = 1 : nspeeds
        for leg = 1 : nlegs
            
            % Get the nuber of samples
            nsamples(trial) = length(sample_cell_list{trial});
            
            % Prealocate
            Fout{ntrials, nspeeds, nlegs} = zeros(size(vars.variables.f, 1), nsamples(trial, speed, leg));
        end
    end
end

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
            for k = sample_cell_list{trial}
                
                % Set model parameters
                model = set_model_parameters(data, vars, model, k, trial, speed, leg);

                % Set model weights
                model = set_model_weights(alpha, vars, model);
                
                % Set cost function normalization
                model = set_model_normalization(data, vars, model);

                % Optimize
%                 fprintf("(%d, %d, %d, %d)\n", k, trial, speed, leg);
                sol = model.solve();

                % Take the solution forces
                f_opt = sol.value(vars.variables.f);

                % Store them in output
                Fout{cnttrials, cntspeeds, cntlegs}(:, cntsamples) = f_opt;
                
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