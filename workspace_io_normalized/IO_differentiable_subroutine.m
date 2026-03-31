function [Fout, Alphaout] = IO_differentiable_subroutine(data, optimal_data, vars, model, sample_list, trial_list, speed_list, leg_list)

% Number of optimizations to do
nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

% Create output structure
n = size(vars.variables.f, 1);      % number of doc vars
ne = size(vars.parameters.A, 1);    % number of eq constr == size of matrix of linear constr
ni = 2*size(vars.variables.f, 1);   % number of ineq constr == nb of bound constr == 2 * nb vars
np = size(vars.variables.alpha);    % number of doc parameters

% Create output structure
Fout = zeros([n, 1, nsamples, ntrials, nspeeds, nlegs]);
Nuout = zeros([ne, 1, nsamples, ntrials, nspeeds, nlegs]);
Lamout = zeros([ni, 1, nsamples, ntrials, nspeeds, nlegs]);
Alphaout = zeros([np, 1, nsamples, ntrials, nspeeds, nlegs]);

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
                model = set_io_model_parameters(data, optimal_data, vars, model, k, trial, speed, leg);
                
                % Set cost function normalization
                model = set_model_normalization(data, vars, model);

                % Optimize
                fprintf("(%d, %d, %d, %d)\n", k, trial, speed, leg);
                sol = model.solve();

                % Take the solution forces
                f_opt = sol.value(vars.variables.f);
                alpha_opt = sol.value(vars.variables.alpha);
                nu_opt = sol.value(vars.variables.nu);
                lam_opt = sol.value(vars.variables.lam);

                % Store them in output
                Fout(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = f_opt;
                Nuout(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = nu_opt;
                Lamout(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = lam_opt;
                Alphaout(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = alpha_opt;
                
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