close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

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
trial_list = 1:10;

TAB = table();

% Do every leg, speed and phase
for leg_list = [1, 2]
for speed_list = [1, 5]
for sample_select = [1, 2]
    % select phase
    if sample_select == 1
        sample_list = 1:4:61;
    else
        sample_list = 61:4:101;
    end
    
    % RMSE and CC vector
    RMSE_vec = zeros(15, 1);
    CC_vec = zeros(15, 1);
    
    % Table
    for ii = 1 : 15

        % Alpha initialize
        alpha = zeros(15, 1);
        alpha(ii) = 1;

        % Compute solution
        Fout_curr = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

        % rmse
        fprintf("The RMSE for L%dS%dP%d with CF%d is : %.4f.\n", leg_list , speed_list, sample_select, ii, rmse(Fout_curr, data.f(:,:,sample_list, trial_list, speed_list, leg_list)));
        
        RMSE_vec(ii) = rmse(Fout_curr, data.f(:,:,sample_list, trial_list, speed_list, leg_list));
        CC_vec(ii) = corr2(reshape(Fout_curr, [], 1), reshape(data.f(:,:,sample_list, trial_list, speed_list, leg_list), [], 1));
        
    end
    TAB.(sprintf('L%dS%dP%d_RMSE', leg_list , speed_list, sample_select)) = RMSE_vec;
    TAB.(sprintf('L%dS%dP%d_CC', leg_list , speed_list, sample_select)) = CC_vec;
end
end
end
writetable(TAB, '..\..\bilevel_optim_results\job_doc_comparison\patient4\table.xlsx');