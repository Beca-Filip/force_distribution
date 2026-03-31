function [fout, sflag, emsg] = DO_single(alpha, vars, model, data, sample, trial, speed, leg)

% Set model parameters
model = set_model_parameters(data, vars, model, sample, trial, speed, leg);
% Set model weights
model = set_model_weights(alpha, vars, model);
% Set cost function normalization
model = set_model_normalization(data, vars, model);

% Set success flag to true
sflag = true;
% Set error message to empty string
emsg = "";

try
    % Optimize
    sol = model.solve();
    % Take the solution forces
    fout = sol.value(vars.variables.f);
    % End
    return
catch optiexception
    % If error, output forces are undefined
    fout = nan(size(model.x));
    % Set flag to false
    sflag = false;
    % Print error identifier and message
    fprintf(1, "Error identifier:\n%s", optiexception.identifier);
    fprintf(1, "Error message:\n%s", optiexception.message);
    % Assign error message
    emsg = optiexception.message;
    % End
    return
end
end