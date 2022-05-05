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

% Create weight vector
alpha_ref = zeros(16, 1);
alpha_ref(2) = 1;

% Calculate Forces
fprintf(strcat('DOC:', repmat('\n', 1, nemptylines)));
Fout = DO_subroutine_normalized(alpha_ref, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% Check VALIDITY of DOC
% Plot results
figure;
hold all;
plot_vector_quantities(sample_list-1, squeeze(data.f(:,:,sample_list,trial_list,speed_list,leg_list)), [], 'LineWidth', 2, 'Color', [0,0,0], 'Marker', 's');
plot_vector_quantities(sample_list-1, squeeze(Fout), [], 'LineWidth', 2, 'Color', [.8,0,0], 'Marker', 'o');
sgtitle('Data vs. Minimum Norm solution.');
legend('Data', 'Opt');

figure;
hold all;
plot(sample_list-1, sum(squeeze(data.f(:,:,sample_list,trial_list,speed_list,leg_list)).^2, 1), 'LineWidth', 2, 'Color', [0,0,0], 'Marker', 's');
plot(sample_list-1, sum(squeeze(Fout).^2, 1), 'LineWidth', 2, 'Color', [.8,0,0], 'Marker', 'o');
title('Norm at each sample of Data vs. Minimum Norm solution.');
legend('Data', 'Opt');

%% Perform IOC

% Set all data to zero, except the DOC result
data.f = zeros(size(data.f));
data.f(:,:,sample_list,trial_list,speed_list,leg_list) = Fout;

% Get inequality activity tolerance
iatol = 1e-6;

% Get IO model
[modelio, varsio] = form_kesh_io_model_normalized(iatol, data, vars, model, sample_list, trial_list, speed_list, leg_list);
modelio.solver('ipopt');

% Solve
fprintf(strcat(repmat('\n', 1, nemptylines), 'IOC:', repmat('\n', 1, nemptylines)));
solio = modelio.solve();

% Extract alpha
alpha_ret = solio.value(varsio.variables.alphaio);

%% Extract the regressor

dPhi = varsio.functions.dJset;
dh = varsio.functions.dh;
% Check dh
if isequal(dh, data.A(:, :, sample_list, trial_list, speed_list, leg_list).')
   fprintf('Constraint gradient check successful.\n') 
end
% Check regressors singular values
R = [dPhi, dh];
[UdPhi, SdPhi, VdPhi] = svd(dPhi);
[Udh, Sdh, Vdh] = svd(dh);
[UR, SR, VR] = svd(R);

[QdPhi, RdPhi] = qr(dPhi);

% Plot
figure;
subplot(3, 1, 1)
stem(log10(diag(SR)), 'LineWidth', 2);
title('SVD of the regressor $[\nabla \Phi \ A^T]$', 'interpreter', 'latex', 'fontsize', 13);
ylabel('$\log_{10}(\sigma_i^2)$', 'interpreter', 'latex', 'fontsize', 15);
xlabel('$i$', 'interpreter', 'latex', 'fontsize', 15);
grid;

subplot(3, 1, 2)
stem(log10(diag(SdPhi)), 'LineWidth', 2);
title('SVD of the CF gradient matrix $\nabla \Phi$', 'interpreter', 'latex', 'fontsize', 13);
ylabel('$\log_{10}(\sigma_i^2)$', 'interpreter', 'latex', 'fontsize', 15);
xlabel('$i$', 'interpreter', 'latex', 'fontsize', 15);
grid;

subplot(3, 1, 3)
stem(log10(diag(Sdh)), 'LineWidth', 2);
title('SVD of the lin. cons matrix $A^T$', 'interpreter', 'latex', 'fontsize', 13);
ylabel('$\log_{10}(\sigma_i^2)$', 'interpreter', 'latex', 'fontsize', 15);
xlabel('$i$', 'interpreter', 'latex', 'fontsize', 15);
grid;

%% 

figure;
stem(log10(abs(diag(RdPhi))), 'LineWidth', 2);
title('QR of the CF gradient matrix $\nabla \Phi$', 'interpreter', 'latex', 'fontsize', 13);
ylabel('$\log_{10}(|r_i|)$', 'interpreter', 'latex', 'fontsize', 15);
xlabel('$i$', 'interpreter', 'latex', 'fontsize', 15);
grid;
%% 

% Compare reference values with retrieved values
barval = [alpha_ref, alpha_ret];
figure;

subplot(3,1,1)
bar(alpha_ref)
title('Reference Weights', 'interpreter', 'latex', 'fontsize', 15);
xticks(1:length(alpha_ref))
xticklabels(num2cell(string(1:length(alpha_ref))))
xlabel(sprintf('$\\phi_i$'), 'interpreter', 'latex', 'fontsize', 15);
ylabel(sprintf('$\\omega_i$'), 'interpreter', 'latex', 'fontsize', 15);

subplot(3,1,2)
bar(alpha_ret,'FaceColor',[1,0,0])
title('Retrieved Weights', 'interpreter', 'latex', 'fontsize', 15);
xticks(1:length(alpha_ref))
xticklabels(num2cell(string(1:length(alpha_ref))))
xlabel(sprintf('$\\phi_i$'), 'interpreter', 'latex', 'fontsize', 15);
ylabel(sprintf('$\\omega_i$'), 'interpreter', 'latex', 'fontsize', 15);

subplot(3,1,3)
bar(VR(1:length(alpha_ref), end),'FaceColor',[1,0,0])
title('SVD Smallest Eigenvector', 'interpreter', 'latex', 'fontsize', 15);
xticks(1:length(alpha_ref))
xticklabels(num2cell(string(1:length(alpha_ref))))
xlabel(sprintf('$\\phi_i$'), 'interpreter', 'latex', 'fontsize', 15);
ylabel(sprintf('$\\omega_i$'), 'interpreter', 'latex', 'fontsize', 15);