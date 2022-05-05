function [modelio, varsio] = form_kesh_io_model_normalized(iatol, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%FORM_KESH_IO_MODEL_NORMALIZED forms the Keshavaraz et al. (2011) inverse
%optimization formulation model using the data.

% Create an optimization model for inverse optimization
modelio = casadi.Opti();

% Problem dimensionality and main variable (number of cost functions in the
% set)
nio = length(vars.functions.Jset);
alphaio = modelio.variable(nio);

% Store main inverse optimization variables
varsio.variables.alphaio = alphaio;

% Determine the number of equality constraints per sample
neqconsio = length(vars.functions.h);


% Get the expressions that will need to be evaluated
% Cost function set gradient matrix
dJsetio = jacobian(vars.functions.Jset, vars.variables.f).';
dhio = jacobian(vars.functions.h, vars.variables.f).';
dgio = jacobian(vars.functions.g, vars.variables.f).';
% Get the total residual
varsio.total_resio = zeros(1, 1, 'like', vars.variables.f);

% Counters for the loop
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
                % Set direct model parameters and normalization
                % Set model parameters
                model = set_model_parameters(data, vars, model, k, trial, speed, leg);
                % Set cost function normalization
                model = set_model_normalization(data, vars, model);
                
                % Create inverse model variables and define functions
                % Create equality mutlipliers
%                 lambdaio = modelio.variable(neqconsio);
%                 varsio.eqmult(cntsamples, cnttrials, cntspeeds, cntlegs).lambda = lambdaio;
                varsio.eqmult(cntsamples, cnttrials, cntspeeds, cntlegs).lambda = modelio.variable(neqconsio);
                lambdaio = varsio.eqmult(cntsamples, cnttrials, cntspeeds, cntlegs).lambda;

                % Determine active inequalities
                activeineqio = determine_active_inequalities_normalized(data, vars, model, k, trial, speed, leg, iatol);

                % Create inequality multipliers
                if ~isempty(activeineqio)
                    % Get number of them
                    nineqconsio = length(activeineqio);
                    % Create and store the variable
                    muio = modelio.variable(nineqconsio);
                    varsio.ineqmult(cntsamples, cnttrials, cntspeeds, cntlegs).mu = muio;
                    % Impose non-negativity constraints on it
                    modelio.subject_to(muio >= 0);
                else                    
                    nineqconsio = 0;
                    varsio.ineqmult(cntsamples, cnttrials, cntspeeds, cntlegs).mu = [];
                end
                
                % Evaluate derivatives at current value
                % Cost function set gradients
                dJset = model.debug.value(dJsetio, {vars.variables.f == data.f(:,:,k,trial,speed,leg)});
                
                if rank(dJset, 1e-6) < size(dJset, 2)
                    % QR
                    [qdJset, rdJset] = qr(dJset);
                    % Detect deficient cf and remove it
                    remind = find(abs(diag(rdJset)) < 1e-6);
                    % Remove those indices
                    dJset(:, remind) = [];
                    alphaio = modelio.variable(nio - length(remind));
                    varsio.variables.alphaio = alphaio;
                    vario.remind = remind;                    
                end
                
                % Equality constraint function gradients
                dh = model.debug.value(dhio, {vars.variables.f == data.f(:,:,k,trial,speed,leg)});
                % Inequality constraint function gradients
                dg = model.debug.value(dgio, {vars.variables.f == data.f(:,:,k,trial,speed,leg)});
                
                % Store the gradients
                varsio.functions(cntsamples, cnttrials, cntspeeds, cntlegs).dJset = dJset;
                varsio.functions(cntsamples, cnttrials, cntspeeds, cntlegs).dh = dh;
                varsio.functions(cntsamples, cnttrials, cntspeeds, cntlegs).dg = dg;
                
                % Calculate the residual differently if inequalities are
                % present
                if ~isempty(activeineqio)
                    % Get the active constraints only
                    dgact = dg(:, activeineqio);
                    % Compute the residual vector
                    gradlagr = dJset * alphaio + dh * lambdaio + dgact * muio;
                else
                    % Compute the residual vector
                    gradlagr = dJset * alphaio + dh * lambdaio;
                end
                % Compute the squared norm of the residual vector
                resio = sum( (gradlagr).^2, 1 );
                % Store the residual
                varsio.functions(cntsamples, cnttrials, cntspeeds, cntlegs).resio = resio;
                % Add to total residual
                varsio.total_resio = varsio.total_resio + resio;
                
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


% Impose constraints on the main inverse optimization variables
modelio.subject_to(alphaio >= 0);
modelio.subject_to(sum(alphaio) == 1);

% Define the cost function to be the total residual
modelio.minimize(varsio.total_resio);
end