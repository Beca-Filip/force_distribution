function [model] = set_io_model_parameters(data, optimal_data, vars, model, k, trial, speed, leg)

% Set the value of the parameters
%         fmin: [35×1 casadi.MX]
%         fmax: [35×1 casadi.MX]
%         pcsa: [35×1 casadi.MX]
%          vmt: [35×1 casadi.MX]
%            M: [35×1 casadi.MX]
%     fpassive: [35×1 casadi.MX]
%           f0: [35×1 casadi.MX]
%            m: [35×1 casadi.MX]
%            r: [35×1 casadi.MX]
%            A: [5×35 casadi.MX]
%            b: [5×1 casadi.MX]

% Target value
model.set_value(vars.parameters.fdata, data.f(:, :, k, trial, speed, leg));

% Time changing parameters
model.set_value(vars.parameters.fmin, data.fmin(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.fmax, data.fmax(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.vmt, data.vmt(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.M, data.M(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.fpassive, data.fpassive(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.A, data.A(:, :, k, trial, speed, leg));
model.set_value(vars.parameters.b, data.b(:, :, k, trial, speed, leg));

% Time constant parameters
model.set_value(vars.parameters.pcsa, data.pcsa(:, :, trial, speed, leg));
model.set_value(vars.parameters.f0, data.f0(:, :, trial, speed, leg));
model.set_value(vars.parameters.m, data.mass(:, :, trial, speed, leg));
model.set_value(vars.parameters.r, data.r(:, :, trial, speed, leg));

% Set initial guess
model.set_initial(vars.variables.f, optimal_data.f(:, :, k, trial, speed, leg));
model.set_initial(vars.variables.nu, optimal_data.nu(:, :, k, trial, speed, leg));
model.set_initial(vars.variables.lam, optimal_data.lam(:, :, k, trial, speed, leg));
model.set_initial(vars.variables.alpha, optimal_data.alpha(:, :, k, trial, speed, leg));
% model.set_initial(vars.variables.f, (data.fmax(:, :, k, trial, speed, leg) + data.fmin(:, :, k, trial, speed, leg))/2);
% model.set_initial(vars.variables.f, data.finit_feas(:,:,k,trial,speed,leg));
end