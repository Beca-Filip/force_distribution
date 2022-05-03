function [model, vars] = form_casadi_model_normalized()

% Create optimization problem instance
model = casadi.Opti();

% Problem dimensionality and variable
n = 35;
f = model.variable(n);

% Number of equality constraints
ne = 5;

% Create model parameters
fmin = model.parameter(n);
fmax = model.parameter(n);
pcsa = model.parameter(n);
vmt = model.parameter(n);
M = model.parameter(n);
fpassive = model.parameter(n);
f0 = model.parameter(n);
m = model.parameter(n);
r = model.parameter(n);
A = model.parameter(ne, n);
b = model.parameter(ne);

% Get all cost functins
Jset = cost_function_set(f, fmin, fmax, pcsa, vmt, M, fpassive, f0, m, r);

% Normalize costs
J_max = model.parameter(length(Jset));
J_min = model.parameter(length(Jset));
for ii = 1 : length(Jset)
    Jset(ii) = (Jset(ii) - J_min(ii)) / (J_max(ii) - J_min(ii));
end

% Get cost function parametrization
alpha = model.parameter(length(Jset));

% Get cost functions and constraints
J = cost_function_normalized(alpha, f, fmin, fmax, pcsa, vmt, M, fpassive, f0, m, r, J_min, J_max);
h = eq_constraint_function(f, A, b);
g = ineq_constraint_function(f, fmin, fmax);

% Add cost and constraints
model.minimize(J);
model.subject_to(h == 0);
model.subject_to(g <= 0);

% Set variables
vars.variables.f = f;
% Set parameters
vars.parameters.alpha = alpha;
vars.parameters.fmin = fmin;
vars.parameters.fmax = fmax;
vars.parameters.pcsa = pcsa;
vars.parameters.vmt = vmt;
vars.parameters.M = M;
vars.parameters.fpassive = fpassive;
vars.parameters.f0 = f0;
vars.parameters.m = m;
vars.parameters.r = r;
vars.parameters.A = A;
vars.parameters.b = b;
vars.parameters.J_min = J_min;
vars.parameters.J_max = J_max;
% Set optimization functions
vars.functions.J = J;
vars.functions.Jset = Jset;
vars.functions.h = h;
vars.functions.g = g;

end