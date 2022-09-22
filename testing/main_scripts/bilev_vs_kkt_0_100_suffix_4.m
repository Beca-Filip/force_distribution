% doc_sample_list = 1:4:61;
doc_sample_list = 61:4:101;
% doc_sample_list = 1:4:101;

% For each cost fun have a Fout structure and a RMSE and CC structure
Fout = cell(1, size(vars.functions.Jset, 2)+2);
RMSE_pmt = cell(1, size(vars.functions.Jset, 2)+2);
TPCC_pmt = cell(1, size(vars.functions.Jset, 2)+2);

% For each cost function
for cf = 1 : size(vars.functions.Jset, 2)
    
    % Curr alpha
    alpha = zeros(size(vars.functions.Jset, 2), 1);
    alpha(cf) = 1;
    % Curr force predictions
    Fcurr = DO_subroutine_normalized(alpha, data, vars, model, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
    Fcurr_aug = data.f;
    Fcurr_aug(:,:,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list) = Fcurr;
    Fout{cf} = Fcurr_aug;
    
    % For each muscle calculate RMSE
    RMSE_pmt{cf} = force_rmse_per_muscle_trajectory(Fout{cf}, data.f,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list);
    % For each muscle calculate TPCC
    TPCC_pmt{cf} = force_correlation_per_muscle_trajectory(Fout{cf}, data.f,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list);
end

%% % Bilevel prediction error
cf = size(vars.functions.Jset, 2) + 1;
Fout{cf} = data_pred_bilev.f;

% For each muscle calculate RMSE
RMSE_pmt{cf} = force_rmse_per_muscle_trajectory(Fout{cf}, data.f,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list);
% For each muscle calculate TPCC
TPCC_pmt{cf} = force_correlation_per_muscle_trajectory(Fout{cf}, data.f,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list);

%% % IKKT prediction error

cf = size(vars.functions.Jset, 2) + 2;
Fout{cf} = data_pred_KKT.f;

% For each muscle calculate RMSE
RMSE_pmt{cf} = force_rmse_per_muscle_trajectory(Fout{cf}, data.f,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list);
% For each muscle calculate TPCC
TPCC_pmt{cf} = force_correlation_per_muscle_trajectory(Fout{cf}, data.f,doc_sample_list,doc_trial_list,doc_speed_list,doc_leg_list);

%% Plot graph

% Mean across trial
mtr_RMSE_pmt = zeros(size(data.f, 1), size(vars.functions.Jset, 2) + 2);
mtr_TPCC_pmt = zeros(size(data.f, 1), size(vars.functions.Jset, 2) + 2);

% Std across trial
str_RMSE_pmt = zeros(size(data.f, 1), size(vars.functions.Jset, 2) + 2);
str_TPCC_pmt = zeros(size(data.f, 1), size(vars.functions.Jset, 2) + 2);

% Std across trials averaged over muscles
mmusc_str_RMSE_pmt = zeros(1, size(vars.functions.Jset, 2) + 2);
mmusc_str_TPCC_pmt = zeros(1, size(vars.functions.Jset, 2) + 2);

% Mean across all
mall_RMSE_pmt = zeros(1, size(vars.functions.Jset, 2) + 2);
mall_TPCC_pmt = zeros(1, size(vars.functions.Jset, 2) + 2);

% Std across all
sall_RMSE_pmt = zeros(1, size(vars.functions.Jset, 2) + 2);
sall_TPCC_pmt = zeros(1, size(vars.functions.Jset, 2) + 2);

% 
for cf = 1 : size(vars.functions.Jset, 2) + 2
    % Calculate per muscle means across trials
    mtr_RMSE_pmt(:, cf) = mean(RMSE_pmt{cf}, 2);
    mtr_TPCC_pmt(:, cf) = mean(TPCC_pmt{cf}, 2);
    % Calculate per muscle stds across trials
    str_RMSE_pmt(:, cf) = std(RMSE_pmt{cf}, 0, 2);
    str_TPCC_pmt(:, cf) = std(TPCC_pmt{cf}, 0, 2);
    
    % Calculate all means across trials and muscles
    mall_RMSE_pmt(cf) = mean(RMSE_pmt{cf}, 'all');
    mall_TPCC_pmt(cf) = mean(TPCC_pmt{cf}, 'all');
    
    % Calculate all stds across trials and muscles
    sall_RMSE_pmt(cf) = std(RMSE_pmt{cf}, 0, 'all');
    sall_TPCC_pmt(cf) = std(TPCC_pmt{cf}, 0, 'all');
end

% std over trials averaged over muscles
mmusc_str_RMSE_pmt = mean(str_RMSE_pmt, 1);
mmusc_str_TPCC_pmt = mean(str_TPCC_pmt, 1);


% Sort by overall mean
% [~, sorti_rmse] = sort(mall_RMSE_pmt, 'ascend');
stance_rmse = [258; 223; 223; 224; 228; 184; 187; 201; 252; 202; 196; 192; 301; 212; 236; 176; 210];
swing_rmse  = [200; 196; 194; 194; 177; 164; 164; 172; 200; 192; 188; 183; 222; 173; 203; 162; 218];
[~, sorti_rmse] = sort(sqrt(stance_rmse.^2+swing_rmse.^2), 'ascend');
[~, sorti_tpcc] = sort(mall_TPCC_pmt, 'descend');

figure;
cmap = linspecer(size(vars.functions.Jset,2)+2);
b_rmse = bar3(mtr_RMSE_pmt(:, sorti_rmse));
for cf = 1 : length(b_rmse)
    b_rmse(cf).FaceColor = cmap(cf, :);
end

xlabel('$\phi$ axis', 'interpreter', 'latex', 'fontsize', 13);
ylabel('Muscle axis', 'interpreter', 'latex', 'fontsize', 13);
zlabel('RMSE axis [N]', 'interpreter', 'latex', 'fontsize', 13);

xtickloc = 1:17;
xticklab = arrayfun(@(ii) sprintf('$\\phi_{%d}$', ii), xtickloc(1:end-2), 'UniformOutput', false);
xticklab{end+1} = '$J_{Bilevel}$';
xticklab{end+1} = '$J_{IKKT}$';
xticklab = xticklab(sorti_rmse);

xticks(xtickloc);
xticklabels(xticklab);
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';

ytickloc = 1:35;
yticklab = arrayfun(@(ii) sprintf('%d', ii), ytickloc, 'UniformOutput', false);

yticks(ytickloc);
yticklabels(yticklab);
yaxisproperties = get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';

xlim([0.5, size(vars.functions.Jset, 2)+2+.5])
ylim([0.5, size(data.f, 1)+0.5])

[az, el] = view;
view(az-30,el+10);

%%


figure;
b_tpcc = bar3(mtr_TPCC_pmt(:, sorti_tpcc));

xlabel('$\phi$ axis', 'interpreter', 'latex', 'fontsize', 13);
ylabel('Muscle axis', 'interpreter', 'latex', 'fontsize', 13);
zlabel('CC axis', 'interpreter', 'latex', 'fontsize', 13);

xtickloc = 1:17;
xticklab = arrayfun(@(ii) sprintf('$\\phi_{%d}$', ii), xtickloc(1:end-2), 'UniformOutput', false);
xticklab{end+1} = '$J_{Bilevel}$';
xticklab{end+1} = '$J_{IKKT}$';
xticklab = xticklab(sorti_rmse);

xticks(xtickloc);
xticklabels(xticklab);
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';

ytickloc = 1:35;
yticklab = arrayfun(@(ii) sprintf('%d', ii), ytickloc, 'UniformOutput', false);

yticks(ytickloc);
yticklabels(yticklab);
yaxisproperties = get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';

xlim([0.5, size(vars.functions.Jset, 2)+2+.5])
ylim([0.5, size(data.f, 1)+0.5])

[az, el] = view;
view(az-30,el+10);

%%

figure;
b_all_rmse = bar([mall_RMSE_pmt(sorti_rmse).', mmusc_str_RMSE_pmt(sorti_rmse).'], 'stacked');

xlabel('$\phi$ axis', 'interpreter', 'latex', 'fontsize', 13);
ylabel('RMSE axis [N]', 'interpreter', 'latex', 'fontsize', 13);

xtickloc = 1:17;
xticklab = arrayfun(@(ii) sprintf('$\\phi_{%d}$', ii), xtickloc(1:end-2), 'UniformOutput', false);
xticklab{end+1} = '$J_{Bilevel}$';
xticklab{end+1} = '$J_{IKKT}$';
xticklab = xticklab(sorti_rmse);

xticks(xtickloc);
xticklabels(xticklab);
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';

xlim([0.5, size(vars.functions.Jset, 2)+2+.5])

%%


figure;
b_all_tpcc = bar([mall_TPCC_pmt(sorti_tpcc).', mmusc_str_TPCC_pmt(sorti_tpcc).'], 'stacked');

xlabel('$\phi$ axis', 'interpreter', 'latex', 'fontsize', 13);
ylabel('CC axis', 'interpreter', 'latex', 'fontsize', 13);

xtickloc = 1:17;
xticklab = arrayfun(@(ii) sprintf('$\\phi_{%d}$', ii), xtickloc(1:end-2), 'UniformOutput', false);
xticklab{end+1} = '$J_{Bilevel}$';
xticklab{end+1} = '$J_{IKKT}$';
xticklab = xticklab(sorti_rmse);

xticks(xtickloc);
xticklabels(xticklab);
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';

xlim([0.5, size(vars.functions.Jset, 2)+1+.5])
