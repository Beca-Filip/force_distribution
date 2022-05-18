function [model, vars] = form_casadi_model(varargin)
%FORM_CASADI_MODEL creates a casadi opti stack with cost functions from the
%basis.
%   
%   [model, vars] = FORM_CASADI_MODEL()
%   [model, vars] = FORM_CASADI_MODEL(cf_exclude)
%
%   Inputs:
%   cf_exclude ~ contains the order of the cost functions to exclude

% Treat inputs
if ~isempty(varargin)
   cf_exclude = varargin{1};
else
    cf_exclude = [];
end

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
Jset = cost_function_set(f, fmin, fmax, pcsa, vmt, M, fpassive, f0, m, r, cf_exclude);

% Get cost function parametrization
alpha = model.parameter(length(Jset));

% Get cost function
J = zeros(1,1,'like',casadi.MX);
for ii = 1 : length(Jset)
    J = J + alpha(ii) * Jset(ii);
end
% Get constraints
h = eq_constraint_function(f, A, b);
g = ineq_constraint_function(f, fmin, fmax);
ceq = h == 0;
c = g <= 0;

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
vars.functions.Jset = Jset;
vars.functions.h = h;
vars.functions.g = g;
vars.functions.ceq = ceq;
vars.functions.c = c;
vars.functions.cf_exclude = cf_exclude;


end