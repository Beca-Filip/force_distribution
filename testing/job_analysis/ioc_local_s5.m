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

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);

% Perform IOC on these
trial_list = 1:10;

TAB = table();

% Do every leg, speed and phase
for leg_list = [1, 2]
for speed_list = [1, 2, 3, 4, 5]
for sample_select = [1, 2]
    if sample_select == 1
        sample_list = 1:61;
    else
        sample_list = 61:101;
    end
    
    alpha0 = ones(15, 1);
    alpha0(7) = 100;
    alpha0 = alpha0 ./ sum(alpha0);

    % Compute solution
    fprintf("Computing IOC %d-%d-%d.\n", sample_select, speed_list, leg_list);
    % Fout_curr = DO_subroutine_normalized(alpha0, data, vars, model, sample_list, trial_list, speed_list, leg_list);
    [alpha, err, exflg, out, lam, grad, hess] = IO_local_fmincon_search_normalized(alpha0, data, vars, model, sample_list, trial_list, speed_list, leg_list);
    
    % Save solution
    save(sprintf("ios5/%d-%d-%d.mat", sample_select, speed_list, leg_list), "alpha", "err", "exflg", "out", "lam", "grad", "hess", "sample_select", "sample_list", "trial_list", "speed_list", "leg_list");
end
end
end