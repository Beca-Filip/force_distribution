function model = set_model_normalization(data, vars, model)

J_min = data.J_min;
J_max = data.J_max;

J_min(vars.functions.cf_exclude) = [];
J_max(vars.functions.cf_exclude) = [];

model.set_value(vars.parameters.J_min, J_min);
model.set_value(vars.parameters.J_max, J_max);
end