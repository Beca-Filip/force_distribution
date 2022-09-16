% Separate RMSE according to individual muscles and individual cost
% functions
RMSE_musc_fun = zeros(size(data.f, 1), size(vars.functions.Jset, 2) + 1, length(doc_trial_list), length(doc_speed_list), length(doc_leg_list));
% Separate TPCC acording to individual muscles and individual cost
% functions
TPCC_musc_fun = zeros(size(data.f, 1), size(vars.functions.Jset, 2) + 1, length(doc_trial_list), length(doc_speed_list), length(doc_leg_list));

% Cost fun parametrization
alpha = zeros(size(vars.functions.Jset, 2), 1);

% For each cost function
for jj = 1 : size(vars.functions.Jset, 2)
    
    % Curr alpha
    alpha = zeros(size(vars.functions.Jset, 2), 1);
    alpha(jj) = 1;
    % Curr force predictions
    Fcurr = DO_subroutine_normalized(alpha, data, vars, model, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list);
    
    % For each muscle calculate RMSE
    for ii = 1 : size(data.f, 1)
        for itr = 1: length(doc_trial_list)
            for isp = 1 : length(doc_speed_list)
                for ile = 1 : length(doc_leg_list)
                    RMSE_musc_fun(ii, jj, itr, isp, ile) = rmse(data.f(ii, :, doc_sample_list, doc_trial_list(itr), doc_speed_list(isp), doc_leg_list(ile)), Fcurr(ii, :, :, itr, isp, ile));
                    TPCC_musc_fun(ii, jj, itr, isp, ile) = corr2(...
                        squeeze(data.f(ii, :, doc_sample_list, doc_trial_list(itr), doc_speed_list(isp), doc_leg_list(ile))), ...
                        squeeze(Fcurr(ii, :, :, itr, isp, ile)));    
                end
            end
        end
    end
end

%% % Bilevel prediction error
jj = size(vars.functions.Jset, 2) + 1;
% For each muscle calculate RMSE
for ii = 1 : size(data.f, 1)
    for itr = 1: length(doc_trial_list)
        for isp = 1 : length(doc_speed_list)
            for ile = 1 : length(doc_leg_list)
                RMSE_musc_fun(ii, jj, itr, isp, ile) = rmse(data.f(ii, :, doc_sample_list, doc_trial_list(itr), doc_speed_list(isp), doc_leg_list(ile)),...
                                                            data_pred_bilev.f(ii, :, doc_sample_list, doc_trial_list(itr), doc_speed_list(isp), doc_leg_list(ile)));
                TPCC_musc_fun(ii, jj, itr, isp, ile) = corr2(squeeze(data.f(ii, :, doc_sample_list, doc_trial_list(itr), doc_speed_list(isp), doc_leg_list(ile))), ...
                                                             squeeze(data_pred_bilev.f(ii, :, doc_sample_list, doc_trial_list(itr), doc_speed_list(isp), doc_leg_list(ile))));
            end
        end
    end
end

%%
all_RMSE_musc_fun = RMSE_musc_fun;
all_TPCC_musc_fun = TPCC_musc_fun;

RMSE_musc_fun = mean(RMSE_musc_fun, [3, 4, 5]);
RMSE_musc_fun_std = std(RMSE_musc_fun, 0, [3, 4, 5]);
TPCC_musc_fun = mean(TPCC_musc_fun, [3, 4, 5]);
TPCC_musc_fun_std = std(TPCC_musc_fun, 0, [3, 4, 5]);

%%
% [~, sorti] = sort(RMSE_musc_fun(33, :));
[~, sorti] = sort(mean(RMSE_musc_fun, 1));

figure;
b = bar3(RMSE_musc_fun(:, sorti));
xlabel('$\phi$ axis', 'interpreter', 'latex');
ylabel('Muscle axis', 'interpreter', 'latex');
zlabel('RMSE axis', 'interpreter', 'latex');


% Tick
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
xaxisproperties.FontSize = 13;
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';
yaxisproperties.FontSize = 13;
zaxisproperties= get(gca, 'ZAxis');
zaxisproperties.FontSize = 13;

xtickloc = 1:16;
xticklab = arrayfun(@(ii) sprintf('$\\phi_{%d}$', ii), xtickloc(1:end-1), 'UniformOutput', false);
xticklab{end+1} = '$J_{\rm Bilevel}$';
xticks(xtickloc);
xticklabels(xticklab(sorti));

ytickloc = 1:35;
yticklab = arrayfun(@(ii) sprintf('%d', ii), ytickloc, 'UniformOutput', false);
yticks(ytickloc);
yticklabels(yticklab);
ylim([0.5, 35.5]);

[az, el] = view;
view(az-30, el+10)

%%


[~, sorti] = sort(mean(TPCC_musc_fun, 1));

figure;
b = bar3(TPCC_musc_fun(:, sorti));
xlabel('$\phi$ axis', 'interpreter', 'latex');
ylabel('Muscle axis', 'interpreter', 'latex');
zlabel('CC axis', 'interpreter', 'latex');


% Tick
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
xaxisproperties.FontSize = 13;
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';
yaxisproperties.FontSize = 13;
zaxisproperties= get(gca, 'ZAxis');
zaxisproperties.FontSize = 13;

xtickloc = 1:16;
xticklab = arrayfun(@(ii) sprintf('$\\phi_{%d}$', ii), xtickloc(1:end-1), 'UniformOutput', false);
xticklab{end+1} = '$J_{\rm Bilevel}$';
xticks(xtickloc);
xticklabels(xticklab(sorti));

ytickloc = 1:35;
yticklab = arrayfun(@(ii) sprintf('%d', ii), ytickloc, 'UniformOutput', false);
yticks(ytickloc);
yticklabels(yticklab);
ylim([0.5, 35.5]);

[az, el] = view;
view(az-30, el+10)

%%

% Mean RMSE
mean_RMSE_musc_fun = mean(RMSE_musc_fun, 1);

% % Concatenate
% conc_RMSE_musc_fun = cat(1, RMSE_musc_fun, zeros(5, size(RMSE_musc_fun, 2)), mean_RMSE_musc_fun);

figure;
b = bar3(mean_RMSE_musc_fun(:, sorti));
ylabel('$\phi$-axis', 'interpreter', 'latex');
zlabel('RMSE-axis', 'interpreter', 'latex');
view(az-90, el);


yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';
yticks(xtickloc);
yticklabels(xticklab(sorti));