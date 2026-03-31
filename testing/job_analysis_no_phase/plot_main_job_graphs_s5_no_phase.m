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

% Which trials to look at 
trial_select = 1 : 10;
sample_select = 1 : 101;

% Do every leg, speed and phase
for leg_select = [1, 2]
for speed_select = 2
    % Load IOC results    
    results_1 = importdata(sprintf('./ios5/whole-cycle-speed-3-leg-%d.mat', leg_select));
    results_2 = importdata(sprintf('./ios4/whole-cycle-speed-2-leg-%d.mat', leg_select));
    
    % Extract alphas
    alpha_1 = results_1.alpha;
    alpha_2 = results_2.alpha;

    % Compute solution
    fprintf("Computing solution speed: %d - leg: %d.\n", speed_select, leg_select);
    % Forces
    predictedForces_1 = DO_subroutine_normalized(alpha_1, data, vars, model, sample_select, trial_select, speed_select, leg_select);
    predictedForces_2 = DO_subroutine_normalized(alpha_2, data, vars, model, sample_select, trial_select, speed_select, leg_select);
    % RMSEs
    predictedRmse_1 = rmse(data.f(:, :, sample_select, trial_select, speed_select, leg_select), predictedForces_1);
    predictedRmse_2 = rmse(data.f(:, :, sample_select, trial_select, speed_select, leg_select), predictedForces_2);
    
    % RMSEs per trial
    predictedRmsePerTrial_1 = nan(size(trial_select));
    predictedRmsePerTrial_2 = nan(size(trial_select));
    for jj = 1 : length(trial_select)
        trial = trial_select(jj);
        predictedRmsePerTrial_1(jj) = rmse(data.f(:, :, sample_select, trial, speed_select, leg_select), predictedForces_1(:, :, sample_select, trial));
        predictedRmsePerTrial_2(jj) = rmse(data.f(:, :, sample_select, trial, speed_select, leg_select), predictedForces_1(:, :, sample_select, trial));
    end
    
    % Plot data
    close all;
    figure('WindowState', 'maximized');
    hold all;
    subjIdx = 2;
    speedIdx = speed_select;
    legIdx = leg_select;
    plotQuantilesAndPredictionsNoPhase(subjIdx, speedIdx, legIdx, data, predictedRmsePerTrial_1, predictedForces_1, predictedForces_2);
    TileFigures;

    exportgraphics(gcf, sprintf('../../bilevel_optim_results/job2_no_phase/subj-%d-force-vs-time-leg-%d.pdf', subjIdx, legIdx), 'ContentType', 'vector');
    exportgraphics(gcf, sprintf('../../bilevel_optim_results/job2_no_phase/subj-%d-force-vs-time-leg-%d.png', subjIdx, legIdx), 'ContentType', 'image', 'Resolution', 300);

end
end