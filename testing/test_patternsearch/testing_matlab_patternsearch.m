% Minimize norm over simplex
clear all;
close all;
clc;

n = 17;
m = 5; % Imagine our gridded simplex with 5 points

% Cost function
fun = @(alpha) 0.5 * alpha.' * alpha;

% Initial solution
alpha0 = zeros(n, 1);
alpha(1) = 1;

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
             'MaxFunctionEvaluations', 1000, ...
             'MaxIterations', 100, ...
             'MaxTime', 10, ...
             'PlotFcn', {'psplotfuncount', 'psplotmaxconstr', 'psplotbestf', 'psplotmeshsize', 'psplotbestx'}, ...
             'InitialMeshSize', mesh_size, ...
             'MeshTolerance', 1e-3, ...
             'UseParallel', true ...
             );
         
% Call pattern search
[alpha_opt, fval_opt, ef_opt, out_opt] = patternsearch(fun, alpha0, A, b, Aeq, beq, lb, ub, nonlcon, ps_options);