function [alpha_opt, fval_opt, varargout] = IO_local_fmincon_search_normalized(alpha0, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%IO_LOCAL_SEARCH performs a local gradient-free search over the cost function
%parametrization.
%
%   [alpha, err] = IO_LOCAL_SEARCH(alpha0, mesh_size, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%   

% Cost function
fun = @(alpha) IO_inner_loop_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% Initial solution
n = length(alpha0);
alpha0 = reshape(alpha0, [], 1);

% Linear constraints
Aeq = ones(1, n);
beq = 1;
A = [];
b = [];

% Bound constraints
lb = zeros(n, 1);
ub = [];

% Nonlinear constraints
nonlcon = [];

% Pattern search options
fmc_options = optimoptions(@fmincon, ...
             'ConstraintTolerance', 1e-6, ...
             'Display', 'iter', ...
             'MaxFunctionEvaluations', 2e4, ...
             'MaxIterations', 1e3, ...
             'PlotFcn', {'optimplotx', 'optimplotfval', 'optimplotconstrviolation', 'optimplotfirstorderopt'}, ...
             'StepTolerance', 1e-6...
         );
         
% Call fmincon
[alpha_opt, fval_opt, ef_opt, out_opt, lambda_opt, grad_opt, hess_opt] = fmincon(fun,alpha0,A,b,Aeq,beq,lb,ub,nonlcon,fmc_options);

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