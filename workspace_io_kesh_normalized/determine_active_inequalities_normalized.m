function act = determine_active_inequalities_normalized(data, vars, model, sample, trial, speed, leg, thresh)
%DETERMINE_ACTIVE_INEQUALITIES_NORMALIZED determines the active inequalities
%at a particular sample,trial,speed,leg quadruple of the data.

% Set model parameters
model = set_model_parameters(data, vars, model, sample, trial, speed, leg);
% Set cost function normalization
model = set_model_normalization(data, vars, model);

% Calculate inequality
g = model.debug.value(vars.functions.g, {vars.variables.f == data.f(:, :, sample, trial, speed, leg)});

% Find active inequalities with desired tolerance
act = find(g >= thresh);

end