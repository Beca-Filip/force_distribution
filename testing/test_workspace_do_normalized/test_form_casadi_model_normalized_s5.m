close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient5.mat';
load(data_dir);

% Add cost function exclusion
cf_exclude = 17;
% Create model
[model, vars] = form_casadi_model_normalized_s5(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);

% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
trial_list = 1:10;
speed_list = 1:5;
leg_list = 1:2;
sample_list = 1:101;

% alpha
alpha = zeros(16, 1);
alpha(1) = 1;

% Start test
vals = [];
for trial = trial_list
for speed = speed_list
for leg = leg_list
for k = sample_list
    model = set_model_parameters(data, vars, model, k, trial, speed, leg);
    model = set_model_weights(alpha, vars, model);
    model = set_model_normalization(data, vars, model);
    
    vals = [vals; ...
        model.debug.value(vars.functions.Jset, {vars.variables.f == data.f(:, :, k, trial, speed, leg)}) ...
        ]; 
end
end
end
end

%%
figure;
hold all;
bar([min(vals); max(vals)].')
