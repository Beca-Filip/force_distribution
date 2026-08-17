function [theta_opt, fval_opt, varargout] = QP_IO_fmincon_search(theta0, data, vars, model, sample_list, trial_list, speed_list, leg_list, cond_max, fmc_overrides)
%QP_IO_FMINCON_SEARCH Local gradient-based IO search over QP cost weights.
%
%   [theta_opt, fval_opt] = QP_IO_FMINCON_SEARCH(theta0, data, vars, model,
%       sample_list, trial_list, speed_list, leg_list, cond_max, fmc_overrides)
%
%   fmc_overrides is an optional struct of fmincon options applied on top of
%   the defaults below, e.g.
%       struct('Display','off', 'PlotFcn',{{}}, 'MaxIterations',50)
%   Required for headless or batch runs: the default PlotFcn opens four
%   figures per fit.
%
%   theta0 is an initial vector of length n*(n+1)/2 + n:
%     theta(1 : n*(n+1)/2)     -> lower-triangular entries of L
%     theta(n*(n+1)/2+1 : end) -> l (n×1 linear weight vector)
%
%   cond_max is optional (default 1e4).  The cost matrix is
%   Q = L*L' + (n/cond_max)*I, and this function adds the equality constraint
%   trace(Q) = n, so that at every FEASIBLE iterate
%
%       cond(Q) <= cond_max
%
%   holds exactly.  Pass Inf to drop both the shift and the constraint, which
%   is the unconstrained formulation -- useful as an ablation.  See
%   THETA_TO_QL for why the constraint has to be an equality (it is what
%   removes the (Q,l) -> (cQ,cl) scale gauge) and QP_TRACE_CONSTRAINT for the
%   residual and its gradient.
%
%   theta0 SHOULD SATISFY the constraint: ||q0||^2 = n*(1 - n/cond_max).
%   L0 = sqrt(1 - n/cond_max) * I gives exactly that, and yields Q0 = I.
%   fmincon tolerates an infeasible start but wastes iterations restoring
%   feasibility.
%
%   Positive-definiteness of Q is guaranteed by Q = L*L' + eps*I.  Sign
%   ambiguity of L is removed by bounding the diagonal entries of L to be
%   non-negative.  No constraints are placed on l.
%
%   To recover Q and l from theta_opt, use the SAME cond_max:
%     [Q, l] = theta_to_Ql(theta_opt, n, cond_max);
%
%   Optional outputs (varargout{1..5}): ef, out, lambda, grad, hess
%   from fmincon.

n = size(vars.variables.f, 1);

if nargin < 9 || isempty(cond_max)
    cond_max = 1e4;
end

if nargin < 10 || isempty(fmc_overrides)
    fmc_overrides = struct();
end

% Cost function
fun = @(theta) QP_IO_inner_loop(theta, data, vars, model, sample_list, trial_list, speed_list, leg_list, cond_max);

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

% Nonlinear equality: trace(Q) = n.  Positive-definiteness itself needs no
% constraint -- it is structural, from Q = L*L' + eps*I.
if isinf(cond_max)
    nonlcon = [];
else
    nonlcon = @(theta) qp_trace_constraint(theta, n, cond_max);
end

% fmincon options
% SpecifyObjectiveGradient:  QP_IO_inner_loop returns [E, dE] when called
%   with two outputs, providing the analytical KKT sensitivity gradient.
% SpecifyConstraintGradient: qp_trace_constraint returns [c, ceq, gc, gceq].
fmc_options = optimoptions(@fmincon, ...
    'SpecifyObjectiveGradient',  true, ...
    'SpecifyConstraintGradient', true, ...
    'ConstraintTolerance',      1e-6, ...
    'Display',                  'iter', ...
    'MaxFunctionEvaluations',   2e4, ...
    'MaxIterations',            1e3, ...
    'PlotFcn',                  {'optimplotx', 'optimplotfval', 'optimplotconstrviolation', 'optimplotfirstorderopt'}, ...
    'StepTolerance',            1e-6 ...
);

% Caller overrides, e.g. struct('Display','off','PlotFcn',{{}},'MaxIterations',50).
% Needed for batch/headless runs -- the PlotFcn above opens four figures per
% fit, which is unusable over 20 fits or on a cluster.
if ~isempty(fmc_overrides)
    override_names = fieldnames(fmc_overrides);
    for oi = 1:numel(override_names)
        fmc_options = optimoptions(fmc_options, ...
            override_names{oi}, fmc_overrides.(override_names{oi}));
    end
end

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
