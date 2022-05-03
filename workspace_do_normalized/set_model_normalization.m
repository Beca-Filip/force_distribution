function model = set_model_normalization(data, vars, model)

model.set_value(vars.parameters.J_min, data.J_min(1:length(vars.parameters.J_min)));
model.set_value(vars.parameters.J_max, data.J_max(1:length(vars.parameters.J_max)));
end