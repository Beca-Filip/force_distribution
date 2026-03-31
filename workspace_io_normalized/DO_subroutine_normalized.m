function [Fout, varargout] = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list)

% Determine if duals should also be computed
if nargout > 1
    dual_flag = true;
else
    dual_flag = false;
end

% Number of optimizations to do
nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

n = size(vars.variables.f, 1);      % number of doc vars
ne = size(vars.parameters.A, 1);    % number of eq constr == size of matrix of linear constr
ni = 2*size(vars.variables.f, 1);   % number of ineq constr == nb of bound constr == 2 * nb vars

% Create output structure
Fout = zeros([n, 1, nsamples, ntrials, nspeeds, nlegs]);
if dual_flag
    Nuout = zeros([ne, 1, nsamples, ntrials, nspeeds, nlegs]);
    Lamout = zeros([ni, 1, nsamples, ntrials, nspeeds, nlegs]);
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
            for k = sample_list

                % Set model parameters
                model = set_model_parameters(data, vars, model, k, trial, speed, leg);

                % Set model weights
                model = set_model_weights(alpha, vars, model);
                
                % Set cost function normalization
                model = set_model_normalization(data, vars, model);

                % Optimize
                % fprintf("(%d, %d, %d, %d)\n", k, trial, speed, leg);
                sol = model.solve();

                % Take the solution forces
                f_opt = sol.value(vars.variables.f);

                if dual_flag
                    dual_vars = sol.value(model.lam_g);
                    Nuout(:, :, cntsamples, trial, speed, leg) = dual_vars(1:ne);
                    Lamout(:, :, cntsamples, trial, speed, leg) = dual_vars(ne+1:ne+ni);
                end

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

if nargout > 1
    varargout{1} = Nuout;
    varargout{2} = Lamout;
end

end