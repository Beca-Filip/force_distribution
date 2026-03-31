function [model, vars] = form_casadi_qp_model(n)
%FORM_CASADI_QP_MODEL Create a CasADi Opti model for a normalized QP.
%   [model, vars] = FORM_CASADI_QP_MODEL(n) constructs a CasADi
%   Opti stack representing a quadratic program expressed in normalized
%   force coordinates. The decision variable is a vector f of length n.
%   The cost is 0.5 * f_norm' * Q * f_norm + l' * f_norm where
%   f_norm = F_invcov * (f - f_mean). Equality and inequality constraints
%   are provided by external functions eq_constraint_function and
%   ineq_constraint_function which accept the (original, not normalized)
%   force f and appropriate parameters.
%
%   Input:
%     n - (optional) dimensionality of decision variable f. Default n = 35.
%
%   The function returns:
%     model - the CasADi Opti object with the problem defined.
%     vars  - a struct collecting handles to variables, parameters and
%             symbolic expressions:
%       vars.variables.f         - decision variable f (n x 1)
%       vars.parameters.*        - parameter handles: f_mean, F_invcov,
%                                  Q, l, fmin, fmax, A, b
%       vars.parameters.f_norm   - normalized force expression
%       vars.functions.*         - symbolic expressions for J, ceq, c
%
%   Notes:
%     - The parameters fmin, fmax, A, b are created as model parameters
%       below and included in vars.parameters for completeness.
%     - eq_constraint_function and ineq_constraint_function must be on the
%       MATLAB path and accept the arguments used below.
%
%   Example:
%     % Create model and inspect returned handles
%     [model, vars] = form_casadi_qp_model();    % uses n = 35
%     [model, vars] = form_casadi_qp_model(10);  % uses n = 10

if nargin < 1 || isempty(n)
    n = 35;
end

% Create optimization problem instance
model = casadi.Opti('conic');

% Problem dimensionality and variable
f = model.variable(n);

% Number of equality constraints
ne = 5;

% Create model parameters
fmin = model.parameter(n);
fmax = model.parameter(n);
f_mean = model.parameter(n);
F_invcov = model.parameter(n, n);
Q = model.parameter(n, n);
l = model.parameter(n);
A = model.parameter(ne, n);
b = model.parameter(ne);

% Normalized forces
f_norm = F_invcov * (f - f_mean);

% Get cost function
J = 0.5 * f_norm.' * Q * f_norm + l.' * f_norm;

% Get constraints
ceq = eq_constraint_function(f, A, b);
c = ineq_constraint_function(f, fmin, fmax);

% Add cost and constraints
model.minimize(J);
model.subject_to(ceq == 0);
model.subject_to(c <= 0);

% Set variables
vars.variables.f = f;
% Set parameters
vars.parameters.f_mean = f_mean;
vars.parameters.F_invcov = F_invcov;
vars.parameters.Q = Q;
vars.parameters.l = l;
vars.parameters.fmin = fmin;
vars.parameters.fmax = fmax;
vars.parameters.A = A;
vars.parameters.b = b;
vars.parameters.f_norm = f_norm;
% Set optimization functions
vars.functions.J = J;
vars.functions.ceq = ceq;
vars.functions.c = c;

end