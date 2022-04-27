function [alpha_opt, fval_opt] = IO_local_search(alpha0, mesh_size, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%IO_LOCAL_SEARCH performs a local gradient-free search over the cost function
%parametrization.
%
%   [alpha, err] = IO_LOCAL_SEARCH(alpha0, mesh_size, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%   

% Cost function
fun = @(alpha) IO_inner_loop(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% Initial solution
n = length(alpha0);

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
ps_options = optimoptions(@patternsearch, ...
             'ConstraintTolerance', 1e-6, ...
             'Display', 'iter', ...
             'MaxFunctionEvaluations', 250, ...
             'MaxIterations', 10, ...
             'PlotFcn', {'psplotfuncount', 'psplotmaxconstr', 'psplotbestf', 'psplotmeshsize', 'psplotbestx'}, ...
             'InitialMeshSize', 1, ...
             'MeshTolerance', 1e-3...
             );
         
% Call pattern search
[alpha_opt, fval_opt, ef_opt, out_opt] = patternsearch(fun, alpha0, A, b, Aeq, beq, lb, ub, nonlcon, ps_options);

end