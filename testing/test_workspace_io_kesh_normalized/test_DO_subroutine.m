close all;
clear all;
clc;

nemptylines = 5;
% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Create DOC model
[model, vars] = form_casadi_model_normalized();
model.solver('ipopt');

% Perform DOC-IOC on these
trial_list = 9;
speed_list = 5;
leg_list = 1;
sample_list = 1;

% Check VALIDITY of DOC
% Plot results
figure(1)
hold all;
plot_vector_quantities(sample_list-1, squeeze(data.f(:,:,sample_list,trial_list,speed_list,leg_list)), [], 'LineWidth', 2, 'Color', [0,0,0], 'Marker', 's', 'DisplayName', 'Data');
sgtitle('Data vs. Optimization solutions.');
legend('interpreter', 'latex', 'location', 'best', 'numcolumns', 5);

% Plot 2-Norm
figure(2)
hold all;
plot(sample_list-1, sum(squeeze(data.f(:,:,sample_list,trial_list,speed_list,leg_list)).^2, 1), 'LineWidth', 2, 'Color', [0,0,0], 'Marker', 's', 'DisplayName', 'Data');
title('Norm at each sample of Data vs. Optimization solutions.');
legend('interpreter', 'latex', 'location', 'best', 'numcolumns', 5);

% Plot Activations
normalize_force = @(f,fmin,fmax) (f-fmin)./(fmax-fmin);
figure(3)
hold all;
plot(1:35, normalize_force(...
                           data.f(:,:,sample_list,trial_list,speed_list,leg_list),...
                           data.fmin(:,:,sample_list,trial_list,speed_list,leg_list),...
                           data.fmax(:,:,sample_list,trial_list,speed_list,leg_list)...
                           ),...
    'LineWidth', 2, 'Color', [0,0,0], 'Marker', 's', 'DisplayName', 'Data');
title('Data Activations vs. Optimization Activations');
legend('interpreter', 'latex', 'location', 'best', 'numcolumns', 5);

% Create a colormap
cmap = linspecer(16);
% Optimize all cost functions
for ii = 1 : 16
    % Reset alpha
    alpha_ref = zeros(16, 1);
    alpha_ref(ii) = 1;

    % Calculate Forces
%     fprintf(strcat('DOC:', repmat('\n', 1, nemptylines)));
    Fout = DO_subroutine_normalized(alpha_ref, data, vars, model, sample_list, trial_list, speed_list, leg_list);
    
    % Plot force vector
    figure(1);
    plot_vector_quantities(sample_list-1, squeeze(Fout), [], 'LineWidth', 2, 'Color', cmap(ii, :), 'Marker', 'o', 'DisplayName', sprintf('$\\phi_{%d}$', ii));
    % Plot force 2 norm
    figure(2);
    plot(sample_list-1, sum(squeeze(Fout).^2, 1), 'LineWidth', 2, 'Color', cmap(ii, :), 'Marker', 'o', 'DisplayName', sprintf('$\\phi_{%d}$', ii));
    % Plot activations
    figure(3)
    plot(1:35, normalize_force(...
                           Fout,...
                           data.fmin(:,:,sample_list,trial_list,speed_list,leg_list),...
                           data.fmax(:,:,sample_list,trial_list,speed_list,leg_list)...
                           ),...
    'LineWidth', 2, 'Color', cmap(ii, :), 'Marker', 'o', 'DisplayName', sprintf('$\\phi_{%d}$', ii));
end

