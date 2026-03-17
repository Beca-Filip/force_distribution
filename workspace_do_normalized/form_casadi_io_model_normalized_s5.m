function [model, vars] = form_casadi_io_model_normalized_s5(varargin)
%FORM_CASADI_IO_MODEL_NORMALIZED_S5 creates a casadi opti stack with cost
%functions from the basis and integrating normalization (specialized for
%subject 5).
%   
%   [model, vars] = FORM_CASADI_IO_MODEL_NORMALIZED_S5()
%   [model, vars] = FORM_CASADI_IO_MODEL_NORMALIZED_S5(cf_exclude)
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
n = 34;
f = model.variable(n);

% Number of equality constraints
ne = 5;
nu = model.variable(ne);
lam = model.variable(2*n);

% Create example
fdata = model.parameter(n);

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

% Normalize costs
J_max = model.parameter(length(Jset));
J_min = model.parameter(length(Jset));
for ii = 1 : length(Jset)
    Jset(ii) = (Jset(ii) - J_min(ii)) / (J_max(ii) - J_min(ii));
end

% Get jacobian
dJset = jacobian(Jset, f);

% Get cost function parametrization
alpha = model.variable(length(Jset));

% Get cost function
J = zeros(1,1,'like',casadi.MX);
for ii = 1 : length(Jset)
    J = J + alpha(ii) * Jset(ii);
end

% Get constraints
h = eq_constraint_function(f, A, b);
g = ineq_constraint_function(f, fmin, fmax);

% Get jacobians
dh = A;
dg = vertcat(-eye(n), eye(n));

% Form KKT
Kstat = dJset.' * alpha + dh.' * nu + dg .' * lam;
Kprim = h;
Kcomp = g .* lam;

% Form constraints
ceq_stat = (Kstat == 0);
ceq_prim = (Kprim == 0);
ceq_comp = (Kcomp == 0);
c_dual = (lam >= 0);
c_prim = (g <= 0);

ceq_norm = (sum(alpha) - 1 == 0);
c_pos = (alpha >= 0);

% Form cost
cost = (f - fdata).' * (f - fdata);

% Add cost and constraints
model.minimize(cost);
model.subject_to(ceq_stat);
model.subject_to(ceq_prim);
model.subject_to(ceq_comp);
model.subject_to(c_dual);
model.subject_to(c_prim);
model.subject_to(ceq_norm);
model.subject_to(c_pos);

% Set variables
vars.variables.f = f;
vars.variables.nu = nu;
vars.variables.lam = lam;
vars.variables.alpha = alpha;
% Set parameters
vars.parameters.fdata = fdata;
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
vars.functions.dJset = dJset;
vars.functions.h = h;
vars.functions.g = g;
vars.functions.ceq_stat = ceq_stat;
vars.functions.ceq_prim = ceq_prim;
vars.functions.ceq_comp = ceq_comp;
vars.functions.c_dual = c_dual;
vars.functions.c_prim = c_prim;
vars.functions.ceq_norm = ceq_norm;
vars.functions.c_pos = c_pos;
vars.functions.cf_exclude = cf_exclude;

end