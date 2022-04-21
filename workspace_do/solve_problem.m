clear all;
close all;
clc;

% Create optimization problem instance
opti = casadi.Opti();

% Problem dimensionality and variable
n = 2;
x = opti.variable(n);

% Call cost function and constraint function
f = cost_function(x);
h = eq_constraint_function(x);
g = ineq_constraint_function(x);

% Add cost and constraints
opti.minimize(f);
opti.subject_to(h == 0);
opti.subject_to(g <= 0);

% Choose solver
opti.solver('ipopt');

% Solve
sol = opti.solve();

x_opt = sol.value(x)
sum(x_opt.^2)