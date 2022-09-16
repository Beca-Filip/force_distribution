% close all
clc

tdpcc = zeros(35, length(doc_sample_list), length(doc_sample_list));

for mm = 1 : 35
    
for itt1 = 1 : length(doc_sample_list)
    for itt2 = 1 : length(doc_sample_list)
        
        tt1 = doc_sample_list(itt1);
        tt2 = doc_sample_list(itt2);
        
        x = data.f(mm, :, tt1, doc_trial_list, doc_speed_list, doc_leg_list);
        y = data_pred_bilev.f(mm, :, tt2, doc_trial_list, doc_speed_list, doc_leg_list);
        x = x(:);
        y = y(:);

        pcc = pearson_correlation_coefficient(x, y);
        
        tdpcc(mm, itt1, itt2) = pcc;

%         figure
%         hold on
%         plot(y);
%         plot(x);
%         xlabel('Realisation');
%         ylabel('Value');
%         legend('y', 'x');
% 
%         figure
%         scatter(x, y, 10, 'x');
    end
end
end

%%

[DOC_SAMPLE_LIST_GRID_1, DOC_SAMPLE_LIST_GRID_2] = meshgrid(doc_sample_list, doc_sample_list);
figrows = 6;
figcols = 6;

fig = figure;
fig.WindowState = 'Maximized';
sgtitle({'Time dependent Pearson cross-correlation coefficient between reference muscle force';
         'trajectories and muscle force trajectory predictions \textbf{for individual muscles}.'},'interpreter','latex');

for ii = 1 : figrows
    for jj = 1 : figcols
        
        curr = (ii-1) * figcols + jj;
        if curr > 35
            break;
        end
        
        subplot(figrows, figcols, curr)
        surf(DOC_SAMPLE_LIST_GRID_1, DOC_SAMPLE_LIST_GRID_2, squeeze(tdpcc(curr, :, :)));
        [az, el] = view;
        hold on;
        
        xlabel('$t_1$', 'interpreter','latex','FontSize',12);
        ylabel('$t_2$', 'interpreter','latex','FontSize',12);
        zlabel(sprintf('$\\rho_{%d}(t_1, t_2)$', curr), 'interpreter','latex','FontSize',12);
%         zlabel(sprintf('$\\rho_{ f_{%d}, \\hat{f}_{%d} }(t_1, t_2)$', curr, curr), 'interpreter','latex');
        plot_x_eq_y_plane([0,100], [-1, 1]);
        view(az, el);
        
        xlim([min(doc_sample_list), max(doc_sample_list)]);
        ylim([min(doc_sample_list), max(doc_sample_list)]);
    end
end

%%

fig = figure;
fig.WindowState = 'Maximized';
sgtitle({'Time dependent Pearson cross-correlation coefficient between reference muscle force';
         'trajectories and muscle force trajectory predictions \textbf{for individual muscles}, for $t_1 = t_2 = t$.'},'interpreter','latex');

for ii = 1 : figrows
    for jj = 1 : figcols
        
        curr = (ii-1) * figcols + jj;
        if curr > 35
            break;
        end
        
        subplot(figrows, figcols, curr)
        plot(doc_sample_list, diag(squeeze(tdpcc(curr, :, :))));
        hold on;
        
        xlabel('$t$', 'interpreter','latex','FontSize',12);
        ylabel('$\rho_{X, Y}(t)$', 'interpreter','latex','FontSize',12);
        
        xlim([min(doc_sample_list), max(doc_sample_list)]); 
    end
end


%% 

mtdpcc = squeeze(mean(tdpcc, 1));

fig = figure;
fig.WindowState = 'Maximized';
sgtitle({'Time dependent Pearson cross-correlation coefficient between reference muscle force';
         'trajectories and muscle force trajectory predictions \textbf{averaged over all muscles}.'},'interpreter','latex');


surf(DOC_SAMPLE_LIST_GRID_1, DOC_SAMPLE_LIST_GRID_2, mtdpcc);
[az, el] = view;
hold on;
plot_x_eq_y_plane([0, 100], [-1, 1]);

xlim([min(doc_sample_list), max(doc_sample_list)]);
ylim([min(doc_sample_list), max(doc_sample_list)]);
view(az, el);
%%

mtdpcc = mean(diag(squeeze(mean(tdpcc, 1))));

fprintf('Mean time dependent Pearson cross-correlation coefficient between reference muscle force\n');
fprintf('trajectories and muscle force trajectory predictions averaged over all muscles, and time,\n');
fprintf('is equal to CC = %.4f.\n', mtdpcc);

fig = figure;
fig.WindowState = 'Maximized';
sgtitle({'Time dependent Pearson cross-correlation coefficient between reference muscle force';
         'trajectories and muscle force trajectory predictions \textbf{averaged over all muscles}, for $t_1 = t_2 = t$.'},'interpreter','latex');


plot(doc_sample_list, diag(squeeze(mean(tdpcc, 1))));
       
xlabel('$t$', 'interpreter','latex','FontSize',12);
ylabel('$\bar{\rho}_{X, Y}(t)$', 'interpreter','latex','FontSize',12);

xlim([min(doc_sample_list), max(doc_sample_list)]);