function model = set_qp_normalization(data, vars, model)

% Set the force normalization parameters (as defined in form_casadi_qp_model.m)
%   f_mean:   [n×1] mean force vector
%   F_invcov: [n×n] inverse square-root covariance (whitening transform)
% The normalized force is: f_norm = F_invcov * (f - f_mean)

model.set_value(vars.parameters.f_mean, data.f_mean);
model.set_value(vars.parameters.F_invcov, data.F_invcov);
end
