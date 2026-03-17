close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient5.mat';
load(data_dir);

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Force bounds
data.fmax(data.fmax <= data.f) = 1.003 * data.f(data.fmax <= data.f);
data.fmin(data.fmin >= data.f) = 0.997 * data.f(data.fmin >= data.f);

epsil = 0.01;
data.fmax(abs(data.fmax - data.f) < epsil) = data.fmax(abs(data.fmax - data.fmin) < epsil) + epsil;
data.fmin(abs(data.f - data.fmin) < epsil) = max(zeros(size(data.fmin(abs(data.f - data.fmin) < epsil))), data.fmin(abs(data.f - data.fmin) < epsil) - epsil);

% Exclude cf
cf_exclude = [16, 17];

% Create model
[model, vars] = form_casadi_model_normalized_s5(cf_exclude);
[modelio, varsio] = form_casadi_io_model_normalized_s5(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);

% Choose DOC init
alpha = zeros(15, 1);
alpha(7) = 1;
sample_list = 1:101;
trial_list = 1:10;
speed_list = 1:5;
leg_list = 1:2;

% Do DOC
[Fout, Nuout, Lamout] = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

%% Create optimal_data struct for IOC

% Choose solver with options
sol_opt= struct;
% sol_opt.ipopt.print_level = 0;
% sol_opt.print_time =0;
% sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
sol_opt.ipopt.max_iter = 2e4;
modelio.solver('ipopt', sol_opt);

optimal_data.f = Fout;
optimal_data.nu = Nuout;
optimal_data.lam = Lamout;
optimal_data.alpha = repmat(alpha, [1, 1, length(sample_list), length(trial_list), length(speed_list), length(leg_list)]);

[Fout, Alphaout, Nuout, Lamout] = IO_subroutine_normalized(data, optimal_data, varsio, modelio, sample_list, trial_list, speed_list, leg_list);
