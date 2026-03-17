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
trial_list = 1:10;
sample_list = 1:101;

TAB = table();

% Do every leg, speed and phase
for leg_list = [1, 2]
for speed_list = 2
    
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
        fprintf("The RMSE for L%dS%d with CF%d is : %.4f.\n", leg_list , speed_list, ii, rmse(Fout_curr, data.f(:,:,sample_list, trial_list, speed_list, leg_list)));
        
        RMSE_vec(ii) = rmse(Fout_curr, data.f(:,:,sample_list, trial_list, speed_list, leg_list));
        CC_vec(ii) = corr2(reshape(Fout_curr, [], 1), reshape(data.f(:,:,sample_list, trial_list, speed_list, leg_list), [], 1));
        
    end
    TAB.(sprintf('L%dS%d_RMSE', leg_list , speed_list)) = RMSE_vec;
    TAB.(sprintf('L%dS%d_CC', leg_list , speed_list)) = CC_vec;
end
end

TAB_IOC = table();

for leg_list = [1, 2]
for speed_list = 2 % must be the same speed list as in the above loop

    % Load IOC results
    s4_ioc_res = load(sprintf('./ios4/whole-cycle-speed-2-leg-%d.mat', leg_list));
    s5_ioc_res = load(sprintf('./ios5/whole-cycle-speed-3-leg-%d.mat', leg_list));

    RMSE_vec = zeros(2, 1);
    CC_vec = zeros(2, 1);
    
    % Perform DOC with S4 ioc result
    alpha = s4_ioc_res.alpha;
    Fout_curr = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);
    fprintf("The RMSE for L%dS%d with the S4 IOC CF is : %.4f.\n", leg_list , speed_list, rmse(Fout_curr, data.f(:,:,sample_list, trial_list, speed_list, leg_list)));
    RMSE_vec(1) = rmse(Fout_curr, data.f(:,:,sample_list, trial_list, speed_list, leg_list));
    CC_vec(1) = corr2(reshape(Fout_curr, [], 1), reshape(data.f(:,:,sample_list, trial_list, speed_list, leg_list), [], 1));

    % Perform DOC with S5 ioc result
    alpha = s5_ioc_res.alpha;
    Fout_curr = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);
    fprintf("The RMSE for L%dS%d with the S5 IOC CF is : %.4f.\n", leg_list , speed_list, rmse(Fout_curr, data.f(:,:,sample_list, trial_list, speed_list, leg_list)));
    RMSE_vec(2) = rmse(Fout_curr, data.f(:,:,sample_list, trial_list, speed_list, leg_list));
    CC_vec(2) = corr2(reshape(Fout_curr, [], 1), reshape(data.f(:,:,sample_list, trial_list, speed_list, leg_list), [], 1));
    
    TAB_IOC.(sprintf('L%dS%d_RMSE', leg_list , speed_list)) = RMSE_vec;
    TAB_IOC.(sprintf('L%dS%d_CC', leg_list , speed_list)) = CC_vec;

end
end

TAB = [TAB; TAB_IOC];

writetable(TAB, '..\..\bilevel_optim_results\job2_no_phase\s4-all-cost-functions.xlsx');

% %%
% 
% TAB_RMSE = TAB(:, 1:2:end);
% TAB_RMSE_L1 = TAB_RMSE(:, 1:end/2);
% TAB_RMSE_L2 = TAB_RMSE(:, end/2+1:end);
% 
% TAB_RMSE_L1_STANCE = TAB_RMSE_L1(:, 1:2:end);
% TAB_RMSE_L1_SWING = TAB_RMSE_L1(:, 2:2:end);
% TAB_RMSE_L2_STANCE = TAB_RMSE_L2(:, 1:2:end);
% TAB_RMSE_L2_SWING = TAB_RMSE_L2(:, 2:2:end);
% 
% fnames = arrayfun(@(n) sprintf("$\\phi_{%d}$", n), 1:15, "UniformOutput", false);
% 
% figure('WindowState','maximized');
% t = tiledlayout(2, 2, "TileSpacing", "tight");
% sgtitle('Subject 1 (i.e. Patient 4)', 'Interpreter', 'latex', 'FontSize', 20);
% 
% ax_ymax = 500;
% nexttile; % subplot(2, 2, 1)
% bar(table2array(TAB_RMSE_L1_STANCE).');
% title('Leg 1, Stance', 'Interpreter', 'latex', 'FontSize', 20);
% % xlabel('Speed (linear from $0.4$ to $0.8 \frac{\rm m}{\rm s}$)', 'Interpreter', 'latex', 'FontSize', 20);
% ylabel('Prediction RMSE [$N$]', 'Interpreter', 'latex', 'FontSize', 20);
% ylim([0, ax_ymax]);
% xaxisproperties = get(gca, 'XAxis'); xaxisproperties.TickLabelInterpreter = "latex"; xaxisproperties.FontSize = 15;
% yaxisproperties = get(gca, 'YAxis'); yaxisproperties.TickLabelInterpreter = "latex"; yaxisproperties.FontSize = 15;
% 
% 
% [~, indexMin] = min(table2array(TAB_RMSE_L1_STANCE), [], 1);
% for j = 1 : size(indexMin, 2)
%     text(j, 1.0*ax_ymax, {sprintf("$\\phi_{%d}$", indexMin(j)); sprintf("%.1f [N]", TAB_RMSE_L1_STANCE{indexMin(j), j})}, 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'Bold', 'HorizontalAlignment','center', 'VerticalAlignment', 'top');
% end
% 
% nexttile; % subplot(2, 2, 2)
% bar(table2array(TAB_RMSE_L2_STANCE).');
% title('Leg 2, Stance', 'Interpreter', 'latex', 'FontSize', 20);
% % xlabel('Speed (linear from $0.4$ to $0.8 \frac{\rm m}{\rm s}$)', 'Interpreter', 'latex', 'FontSize', 20);
% % ylabel('Prediction RMSE [$N$]', 'Interpreter', 'latex', 'FontSize', 20);
% ylim([0, ax_ymax]);
% xaxisproperties = get(gca, 'XAxis'); xaxisproperties.TickLabelInterpreter = "latex"; xaxisproperties.FontSize = 15;
% yaxisproperties = get(gca, 'YAxis'); yaxisproperties.TickLabelInterpreter = "latex"; yaxisproperties.FontSize = 15;
% 
% [~, indexMin] = min(table2array(TAB_RMSE_L2_STANCE), [], 1);
% for j = 1 : size(indexMin, 2)
%     text(j, 1.0*ax_ymax, {sprintf("$\\phi_{%d}$", indexMin(j)); sprintf("%.1f [N]", TAB_RMSE_L2_STANCE{indexMin(j), j})}, 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'Bold', 'HorizontalAlignment','center', 'VerticalAlignment', 'top');
% end
% 
% nexttile; % subplot(2, 2, 3)
% bar(table2array(TAB_RMSE_L1_SWING).');
% title('Leg 1, Swing', 'Interpreter', 'latex', 'FontSize', 20);
% xlabel('Speed (linear from $0.4$ to $0.8 \frac{\rm m}{\rm s}$)', 'Interpreter', 'latex', 'FontSize', 20);
% ylabel('Prediction RMSE [$N$]', 'Interpreter', 'latex', 'FontSize', 20);
% ylim([0, ax_ymax]);
% xaxisproperties = get(gca, 'XAxis'); xaxisproperties.TickLabelInterpreter = "latex"; xaxisproperties.FontSize = 15;
% yaxisproperties = get(gca, 'YAxis'); yaxisproperties.TickLabelInterpreter = "latex"; yaxisproperties.FontSize = 15;
% 
% [~, indexMin] = min(table2array(TAB_RMSE_L1_SWING), [], 1);
% for j = 1 : size(indexMin, 2)
%     text(j, 1.0*ax_ymax, {sprintf("$\\phi_{%d}$", indexMin(j)); sprintf("%.1f [N]", TAB_RMSE_L1_SWING{indexMin(j), j})}, 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'Bold', 'HorizontalAlignment','center', 'VerticalAlignment', 'top');
% end
% 
% nexttile; % subplot(2, 2, 4)
% bar(table2array(TAB_RMSE_L2_SWING).');
% title('Leg 2, Swing', 'Interpreter', 'latex', 'FontSize', 20);
% xlabel('Speed (linear from $0.4$ to $0.8 \frac{\rm m}{\rm s}$)', 'Interpreter', 'latex', 'FontSize', 20);
% % ylabel('Prediction RMSE [$N$]', 'Interpreter', 'latex', 'FontSize', 20);
% legend(fnames, 'Interpreter', 'latex', 'FontSize', 15, 'Position', [0.9143 0.2737 0.0541 0.5361]);
% ylim([0, ax_ymax]);
% xaxisproperties = get(gca, 'XAxis'); xaxisproperties.TickLabelInterpreter = "latex"; xaxisproperties.FontSize = 15;
% yaxisproperties = get(gca, 'YAxis'); yaxisproperties.TickLabelInterpreter = "latex"; yaxisproperties.FontSize = 15;
% 
% [~, indexMin] = min(table2array(TAB_RMSE_L2_SWING), [], 1);
% for j = 1 : size(indexMin, 2)
%     text(j, 1.0*ax_ymax, {sprintf("$\\phi_{%d}$", indexMin(j)); sprintf("%.1f [N]", TAB_RMSE_L2_SWING{indexMin(j), j})}, 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'Bold', 'HorizontalAlignment','center', 'VerticalAlignment', 'top');
% end
% 
% exportgraphics(gcf, '..\..\bilevel_optim_results\job_doc_comparison\patient4\table_graphic.pdf', 'ContentType', 'vector');