function [model] = set_qp_weights(Q, l, vars, model)

% Set the QP cost weights (as defined in form_casadi_qp_model.m)
%   Q: [n×n] quadratic weight matrix
%   l: [n×1] linear weight vector
% The cost is: 0.5 * f_norm' * Q * f_norm + l' * f_norm

model.set_value(vars.parameters.Q, Q);
model.set_value(vars.parameters.l, l);
end
