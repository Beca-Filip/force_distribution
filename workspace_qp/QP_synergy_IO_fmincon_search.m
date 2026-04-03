function [theta_opt, fval_opt, varargout] = QP_synergy_IO_fmincon_search(theta0, syn_info, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%QP_SYNERGY_IO_FMINCON_SEARCH  Local gradient-based IO search over synergy
%                               constrained QP cost weights.
%
%   [theta_opt, fval_opt] = QP_SYNERGY_IO_FMINCON_SEARCH(theta0, syn_info, ...)
%
%   theta0 is the initial synergy parameter vector of length syn_info.n_theta:
%     theta0(1 : n_theta_q)     -> free block-diagonal L entries (Q = L*L')
%     theta0(n_theta_q+1 : end) -> l (linear weight vector)
%
%   To recover Q and l from theta_opt, call:
%     [Q, L, ~, l] = synergy_theta_to_Q(theta_opt, syn_info);
%
%   Optional outputs (varargout{1..5}): ef, out, lambda, grad, hess.

fun = @(th) QP_synergy_IO_inner_loop(th, syn_info, data, vars, model, ...
    sample_list, trial_list, speed_list, leg_list);

theta0 = reshape(theta0, [], 1);

% Bound constraints: diagonal entries of each block L must be >= 0
% (removes sign ambiguity); all other theta entries are unbounded.
lb           = -inf(syn_info.n_theta, 1);
lb(syn_info.diag_mask) = 0;
ub = [];

fmc_options = optimoptions(@fmincon, ...
    'SpecifyObjectiveGradient', true, ...
    'ConstraintTolerance',      1e-6, ...
    'Display',                  'iter', ...
    'MaxFunctionEvaluations',   2e4, ...
    'MaxIterations',            1e3, ...
    'PlotFcn',                  {'optimplotx', 'optimplotfval', 'optimplotconstrviolation', 'optimplotfirstorderopt'}, ...
    'StepTolerance',            1e-6 ...
);

[theta_opt, fval_opt, ef_opt, out_opt, lambda_opt, grad_opt, hess_opt] = ...
    fmincon(fun, theta0, [], [], [], [], lb, ub, [], fmc_options);

if nargout > 6, varargout{5} = hess_opt;  end
if nargout > 5, varargout{4} = grad_opt;  end
if nargout > 4, varargout{3} = lambda_opt; end
if nargout > 3, varargout{2} = out_opt;   end
if nargout > 2, varargout{1} = ef_opt;    end

end
