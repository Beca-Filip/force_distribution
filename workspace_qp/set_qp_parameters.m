function [model] = set_qp_parameters(data, vars, model, k, trial, speed, leg)

% Set the value of the QP parameters (as defined in form_casadi_qp_model.m)
%     fmin: [n×1 casadi.MX]  - lower force bounds
%     fmax: [n×1 casadi.MX]  - upper force bounds
%        A: [5×n casadi.MX]  - equality constraint matrix
%        b: [5×1 casadi.MX]  - equality constraint rhs

% Time changing parameters
model.set_value(vars.parameters.fmin, data.fmin(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.fmax, data.fmax(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.A, data.A(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.b, data.b(:, :, k, trial, speed, leg));

% Set initial guess
model.set_initial(vars.variables.f, data.f(:, :, k, trial, speed, leg));
end
