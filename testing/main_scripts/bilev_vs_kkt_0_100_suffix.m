% Separate RMSE according to individual muscles and individual cost
% functions
RMSE_musc_fun = zeros(size(data.f, 1), size(vars.functions.Jset, 2) + 1);

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
        RMSE_musc_fun(ii, jj) = rmse(data.f(ii, :, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list), Fcurr(ii, :, :, :, :, :));
    end
end

% Bilevel prediction error
jj = size(vars.functions.Jset, 2) + 1;
% For each muscle calculate RMSE
for ii = 1 : size(data.f, 1)
    RMSE_musc_fun(ii, jj) = rmse(data.f(ii, :, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list),...
                                 data_pred_bilev.f(ii, :, doc_sample_list, doc_trial_list, doc_speed_list, doc_leg_list));
end

%%

[~, sorti] = sort(RMSE_musc_fun(33, :));

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