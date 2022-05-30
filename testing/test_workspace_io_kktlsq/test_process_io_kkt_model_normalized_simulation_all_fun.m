close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient4.mat';
load(data_dir);

% Exclude cost function 17
cf_exclude = 17;

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);
model.solver('ipopt');

% Perform IOC on these
trial_list = 1:10;
speed_list = 5;
leg_list = 1;
sample_list = 1:4:61;

% Cost function parametrization
% alpha = normalize_vector(rand(16, 1), 1);
alpha = zeros(16, 1);
alpha(16) = 1;

% Perform the DO on these
Fout = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);
% Store the DO in the data
data.f(:, :, sample_list, trial_list, speed_list, leg_list) = Fout;

% Active inequality thresh
aithresh = -1e-5;

% Create IO model
modelio = form_io_kktlsq_struct_normalized(data, model, vars, sample_list, trial_list, speed_list, leg_list);

% Process IO model
modelio = process_io_kkt_model_normalized(modelio, -1e-1);

% Solve using LSQ
solio = solve_io_kkt_model_normalized(modelio);

% Analyze IO model
[~, ic_cf] = indexing_function(modelio, 'cf_parameters', [], [], [], []);
min(svd(modelio.regressor(:, ic_cf)))
min(svd(modelio.regressor))
[u,s,v] = svd(modelio.regressor);


solv.z = v(:, end);
alpha.'
[~, solv.alpha] = positive_negative_orthant_check(normalize_vector(v(1:length(alpha), end), 1));
solv.alpha.'
solio.alpha.'

