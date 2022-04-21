function [model] = set_model_weights(alpha, vars, model)

% Set the value of the function weights

model.set_value(vars.parameters.alpha, alpha);
end