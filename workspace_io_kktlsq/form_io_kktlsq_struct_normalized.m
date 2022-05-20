function modelio = form_io_kktlsq_struct_normalized(data, model, vars, sample_list, trial_list, speed_list, leg_list)
%FORM_IO_KKTLSQ_STRUCT_NORMALIZED forms the Keshavaraz et al. (2011) inverse
%optimization formulation model using the data, stores it inside a struct
%which contains numerical matrices needed for an LSQ resolution.
%
%   modelio = FORM_IO_KKTLSQ_STRUCT_NORMALIZED(data, model, vars, sample_list, trial_list, speed_list, leg_list)


% Create a structure to store data
modelio = struct();

% Put the DOC model and variables as fields
modelio.do_model = model;
modelio.do_vars = vars;

% Number of DOC variables
n = size(vars.variables.f, 1);
modelio.do_n = n;

% Store the lists
modelio.sample_list = sample_list;
modelio.trial_list = trial_list;
modelio.speed_list = speed_list;
modelio.leg_list = leg_list;

% Store number of emements
modelio.num_samples = numel(sample_list);
modelio.num_trials = numel(trial_list);
modelio.num_speeds = numel(speed_list);
modelio.num_legs = numel(leg_list);

% Store total number of data points
modelio.num_total_data_pts = numel(sample_list) * numel(trial_list) * numel(speed_list) * numel(leg_list);

% Problem dimensionality and main variable (number of cost functions in the
% set)
ncf = length(vars.functions.Jset);
modelio.num_cf = ncf;

% Equality constraint dimensionality
nec = length(vars.functions.h);
modelio.num_ec = nec;

% Inequality constraint dimensionality
nic = length(vars.functions.g);
modelio.num_ic = nic;

% Get the expressions that will need to be evaluated
% Cost function set gradient matrix
dJsetio = jacobian(vars.functions.Jset, vars.variables.f).';
dhio = jacobian(vars.functions.h, vars.variables.f).';
dgio = jacobian(vars.functions.g, vars.variables.f).';

% Counters for the loop
cntsamples = 1;
cnttrials = 1;
cntspeeds = 1;
cntlegs = 1;

% Prealocate
% modelio.aic = cell(numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));
% modelio.num_aic = zeros(numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));
modelio.gradients_cf = zeros(n, ncf, numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));
modelio.gradients_ec = zeros(n, nec, numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));
modelio.gradients_ic = zeros(n, nic, numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));
modelio.values_ic = zeros(nic, numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));


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
                
                % Evaluate inequalities at current value
                g = model.debug.value(vars.functions.g, {vars.variables.f == data.f(:,:,k,trial,speed,leg)});
                % Evaluate derivatives at current value
                % Cost function set gradients
                dJset = model.debug.value(dJsetio, {vars.variables.f == data.f(:,:,k,trial,speed,leg)});
                % Equality constraint function gradients
                dh = model.debug.value(dhio, {vars.variables.f == data.f(:,:,k,trial,speed,leg)});
                % Inequality constraint function gradients
                dg = model.debug.value(dgio, {vars.variables.f == data.f(:,:,k,trial,speed,leg)});
                
%                 % Determine active inequality constraints and store them
%                 aic = determine_active_inequalities_normalized(data, vars, model, k, trial, speed, leg, aithresh);
%                 modelio.aic{cntsamples, cnttrials, cntspeeds, cntlegs} = aic;
% 
%                 % Get number of active inequality constraints
%                 naic = length(aic);
%                 modelio.num_aic(cntsamples, cnttrials, cntspeeds, cntlegs) = naic;
%                 
%                 % Extract active inequality gradients
%                 dgaic = dg(:, aic);

                % Store the numerical gradient matrices
                modelio.gradients_cf(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = dJset;
                modelio.gradients_ec(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = dh;
                modelio.gradients_ic(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = dg;
                
                % Store the inequality constraints
                modelio.values_ic(:, cntsamples, cnttrials, cntspeeds, cntlegs) = g;
                
%                 % Calculate current matrix size
%                 s1_dJset = size(dJset, 1);
%                 [s1_grad, s2_grad] = size(modelio.gradients_whole);
%                 
%                 % Augment the matrix
%                 modelio.gradients_whole = [[matrix.gradients_whole, zeros(s1_grad, nec + naic)];
%                                            [dJset, zeros(s1_dJset, s2_grad - ncf), dh, dgaic]];

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