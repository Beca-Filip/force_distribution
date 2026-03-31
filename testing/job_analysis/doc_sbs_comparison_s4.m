close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
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
[model, vars] = form_casadi_model_normalized(cf_exclude);

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
sample_list = 1:101;
trial_list = 1:10;
speed_list = 1:5;
leg_list = 1:2;

% Output prealocation
ncf = 15;
sz_f = size(data.f);
Fout = zeros([sz_f, ncf]);
MSEout = zeros([sz_f(2:end), ncf]);

% Table
for ii = 1 : ncf

    % Alpha initialize
    alpha = zeros(ncf, 1);
    alpha(ii) = 1;

    % Compute solution
    Fout_curr = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);
    
    % Compute MSE
    MSEout_curr = squeeze(sum((Fout_curr - data.f).^2, 1));

    % Store solution
    Fout(:, :, :, :, :, :, ii) = Fout_curr;

    % Store MSE
    MSEout(:, :, :, :, :, ii) = MSEout_curr;


    % rmse
    fprintf("The RMSE with CF%d is : %.4f.\n", ii, sqrt(sum(MSEout ./ numel(Fout_curr), 'all')));
end


save('..\..\bilevel_optim_results\job_doc_comparison\patient4\doc_sbs_comparison_s4--variables.mat');

%%

leg = 2;
speed = 5;
trial = 1;

for ii = 1 : ncf

    figure('WindowState','maximized')
    hold all;
    plot_vector_quantities_opts_shape(0:100, data.fmax(:, :, :, trial, speed, leg), [], [], [6, 6], 'Color', [0.1, 0.1, 0.1]);
    plot_vector_quantities_opts_shape(0:100, data.fmin(:, :, :, trial, speed, leg), [], [], [6, 6], 'Color', [0., 0., 0.]);
    plot_vector_quantities_opts_shape(0:100, data.f(:, :, :, trial, speed, leg), [], [], [6, 6], 'Color', [0., 0., 1.], 'LineWidth', 2);
    plot_vector_quantities_opts_shape(0:100, Fout(:, :, :, trial, speed, leg, ii), [], [], [6, 6], '--', 'Color', [1., 0., 0.], 'LineWidth', 2);
end