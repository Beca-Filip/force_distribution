function [model, vars] = form_casadi_model()

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

% Get cost function parametrization
alpha = model.parameter(length(Jset));

% Get cost functions and constraints
J = cost_function(alpha, f, fmin, fmax, pcsa, vmt, M, fpassive, f0, m, r);
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
% Set optimization functions
vars.functions.J = J;
vars.functions.Jset = Jset;
vars.functions.h = h;
vars.functions.g = g;

end