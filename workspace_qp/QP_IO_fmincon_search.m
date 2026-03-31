function [theta_opt, fval_opt, varargout] = QP_IO_fmincon_search(theta0, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%QP_IO_FMINCON_SEARCH Local gradient-based IO search over QP cost weights.
%
%   [theta_opt, fval_opt] = QP_IO_FMINCON_SEARCH(theta0, data, vars, model,
%       sample_list, trial_list, speed_list, leg_list)
%
%   theta0 is an initial vector of length n*(n+1)/2 + n:
%     theta(1 : n*(n+1)/2)     -> lower-triangular entries of L, where Q = L*L'
%     theta(n*(n+1)/2+1 : end) -> l (n×1 linear weight vector)
%
%   Positive-semidefiniteness of Q is guaranteed by the Cholesky
%   parameterisation Q = L*L'.  Sign ambiguity of L is removed by bounding
%   the diagonal entries of L to be non-negative.  No constraints are placed
%   on l.
%
%   To recover Q and l from theta_opt, call:
%     [Q, L] = chol_vec_to_Q(theta_opt(1:n*(n+1)/2), n);
%     l      = theta_opt(n*(n+1)/2+1:end);
%
%   Optional outputs (varargout{1..5}): ef, out, lambda, grad, hess
%   from fmincon.

n = size(vars.variables.f, 1);

% Cost function
fun = @(theta) QP_IO_inner_loop(theta, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% Initial solution
theta0 = reshape(theta0, [], 1);

% Linear constraints: none
A   = [];
b   = [];
Aeq = [];
beq = [];

% Bound constraints:
%   diagonal entries of L (positions in the lower-tri vector where row == col)
%   are bounded below by 0 to remove the sign ambiguity of L;
%   all other entries of theta are unbounded.
mask = tril(true(n));
[r, c] = find(mask);                    % row/col of each lower-tri entry in column-major order
diag_idx = find(r == c);               % positions in q corresponding to L(k,k)

lb = -inf(length(theta0), 1);
lb(diag_idx) = 0;
ub = [];

% No nonlinear constraints: PSD is guaranteed by Q = L*L'
nonlcon = [];

% fmincon options
fmc_options = optimoptions(@fmincon, ...
    'ConstraintTolerance',    1e-6, ...
    'Display',                'iter', ...
    'MaxFunctionEvaluations', 2e4, ...
    'MaxIterations',          1e3, ...
    'PlotFcn',                {'optimplotx', 'optimplotfval', 'optimplotconstrviolation', 'optimplotfirstorderopt'}, ...
    'StepTolerance',          1e-6 ...
);

% Call fmincon
[theta_opt, fval_opt, ef_opt, out_opt, lambda_opt, grad_opt, hess_opt] = ...
    fmincon(fun, theta0, A, b, Aeq, beq, lb, ub, nonlcon, fmc_options);

if nargout > 6
    varargout{5} = hess_opt;
end
if nargout > 5
    varargout{4} = grad_opt;
end
if nargout > 4
    varargout{3} = lambda_opt;
end
if nargout > 3
    varargout{2} = out_opt;
end
if nargout > 2
    varargout{1} = ef_opt;
end

end
