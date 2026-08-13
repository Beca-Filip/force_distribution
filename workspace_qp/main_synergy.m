close all; clear all; clc;

% =========================================================================
% QP synergy-constrained inverse-optimal-control  --  local fmincon search
%
% Restricts Q to be block-diagonal where each block corresponds to a muscle
% synergy group.  Q_{ij} is free only when muscles i and j belong to the
% same synergy; all cross-synergy entries are identically zero by construction
% (block-diagonal Cholesky: Q = blkdiag(L_1*L_1', ..., L_K*L_K', diag(.))).
%
% The synergy groups are defined by name patterns below and are resolved
% to per-patient muscle indices via build_synergy_info.
%
% Results are saved to results_qp_synergy/ as:
%   subject-{S}-speed-{speed}-leg-{leg}.mat  (weights + metadata)
%   subject-{S}-speed-{speed}-leg-{leg}.png  (prediction comparison)
%   subject-{S}-speed-{speed}-leg-{leg}-weights.png
%   subject-{S}-speed-{speed}-leg-{leg}-eig.png
% =========================================================================

% ---- Synergy group definitions (edit here to change grouping) -----------
% Each cell contains partial name patterns (case-insensitive).
% A muscle is assigned to the first group whose pattern matches its name.
% NOTE: recfem is listed in Synergy 3 (hip flexors, full role); it is NOT
% in Synergy 1, where it has only a partial role per the anatomical literature.
synergy_groups = { ...
    {'glmax1','glmax2','glmax3','vasmed','vasint','vaslat'}, ...  % S1: hip+knee ext
    {'soleus','gasmed','gaslat'}, ...                             % S2: plantarflexors
    {'iliacus','psoas','recfem'}, ...                             % S3: hip flexors
    {'semimem','semiten','bflh','bfsh'}, ...                      % S4: hamstrings
    {'tibant','edl'}, ...                                         % S5: dorsiflexors
    {'glmed1','glmed2','glmed3','glmin1','glmin2','glmin3', ...   % S6: frontal plane
     'addbrev','addlong','addmagDist','addmagIsch','addmagMid','addmagProx'}, ...
};

% ---- Subject definitions ------------------------------------------------
subjects(1).id           = 4;
subjects(1).data_dir     = fullfile('..', 'Optimization Model Data', 'Patient4.mat');
subjects(1).muscle_names_file = fullfile('..', 'patient_4_muscle_names.mat');
subjects(1).speed_list   = 1:5;

subjects(2).id           = 5;
subjects(2).data_dir     = fullfile('..', 'Optimization Model Data', 'Patient5.mat');
subjects(2).muscle_names_file = fullfile('..', 'patient_5_muscle_names.mat');
subjects(2).speed_list   = 1:5;

% ---- Output directory ---------------------------------------------------
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results_qp_synergy');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% =========================================================================
% Main loop over subjects
% =========================================================================
for si = 1:length(subjects)

    subject_id = subjects(si).id;
    fprintf('=== Subject %d ===\n', subject_id);

    % ---- Load data and muscle names -------------------------------------
    load(subjects(si).data_dir);

    mn_struct    = load(subjects(si).muscle_names_file);
    mn_field     = fieldnames(mn_struct);
    muscle_names = mn_struct.(mn_field{1});   % cell array of char vectors

    % ---- Fix force bounds -----------------------------------------------
    data.fmax(data.fmax <= data.f) = 1.003 * data.f(data.fmax <= data.f);
    data.fmin(data.fmin >= data.f) = 0.997 * data.f(data.fmin >= data.f);

    epsil = 0.01;
    data.fmax(abs(data.fmax - data.f) < epsil) = ...
        data.fmax(abs(data.fmax - data.fmin) < epsil) + epsil;
    data.fmin(abs(data.f - data.fmin) < epsil) = max( ...
        zeros(size(data.fmin(abs(data.f - data.fmin) < epsil))), ...
        data.fmin(abs(data.f - data.fmin) < epsil) - epsil);

    % ---- Force-space normalisation --------------------------------------
    data.f_mean   = compute_f_mean(data.f);
    data.F_invcov = compute_F_invcov(data.f, data.f_mean);

    % ---- Build synergy info (patient-specific muscle ordering) ----------
    n        = size(data.f, 1);
    syn_info = build_synergy_info(muscle_names, synergy_groups);

    % Sanity-check: syn_info must match the data dimensionality
    if syn_info.n ~= n
        error('Muscle count mismatch: data has %d muscles, names file has %d.', n, syn_info.n);
    end

    fprintf('  Synergy parameterisation: %d free Q params + %d l params = %d total\n', ...
        syn_info.n_theta_q, n, syn_info.n_theta);
    for k = 1:syn_info.n_groups
        gi = syn_info.group_indices{k};
        fprintf('    Group %d (%d muscles): %s\n', k, numel(gi), ...
            strjoin(muscle_names(gi), ', '));
    end
    si_idx = syn_info.singleton_indices;
    if ~isempty(si_idx)
        fprintf('    Singletons (%d muscles): %s\n', numel(si_idx), ...
            strjoin(muscle_names(si_idx), ', '));
    end

    % ---- Build QP model -------------------------------------------------
    [model, vars] = form_casadi_qp_model(n);
    sol_opt = struct;
    sol_opt.osqp.verbose = 0;
    model.solver('osqp', sol_opt);

    % ---- Lists ----------------------------------------------------------
    trial_list  = 1:10;
    sample_list = 1:101;

    % ---- Initial theta_syn: block L = identity-like, l = random ---------
    % For each group block of size m_k: diagonal entries of L_k = 1,
    % off-diagonal entries = 0.  Singletons: entry = 1.
    % This yields Q = I restricted to the synergy pattern.
    q_syn0 = zeros(syn_info.n_theta_q, 1);
    q_syn0(syn_info.diag_mask) = 1;   % diagonal entries -> 1
    l0     = randn(n, 1);
    theta0 = [q_syn0; l0];

    % =========================================================================
    % Loop over speeds and legs
    % =========================================================================
    for speed_list = subjects(si).speed_list
    for leg_list = [1, 2]

        fprintf('  Subject %d  Speed %d  Leg %d\n', subject_id, speed_list, leg_list);

        % ---- IO search --------------------------------------------------
        [theta_opt, fval_opt, ef_opt, out_opt, lambda_opt, grad_opt, hess_opt] = ...
            QP_synergy_IO_fmincon_search(theta0, syn_info, data, vars, model, ...
                sample_list, trial_list, speed_list, leg_list);

        % ---- Reconstruct Q and l ----------------------------------------
        [Q_opt, L_opt, ~, l_opt] = synergy_theta_to_Q(theta_opt, syn_info);

        % ---- Predicted forces -------------------------------------------
        Fout = QP_subroutine(Q_opt, l_opt, data, vars, model, ...
            sample_list, trial_list, speed_list, leg_list);

        % ---- Per-trial RMSE ---------------------------------------------
        ntrials = length(trial_list);
        rmse_per_trial = zeros(ntrials, 1);
        for ti = 1:ntrials
            f_ref_t  = data.f(:, :, sample_list, trial_list(ti), speed_list, leg_list);
            f_pred_t = Fout(:, :, :, ti, 1, 1);
            rmse_per_trial(ti) = rmse(f_ref_t(:), f_pred_t(:));
        end

        % ---- Save .mat --------------------------------------------------
        fname_base = fullfile(results_dir, ...
            sprintf('subject-%d-speed-%d-leg-%d', subject_id, speed_list, leg_list));

        save([fname_base '.mat'], ...
            'theta_opt', 'Q_opt', 'L_opt', 'l_opt', 'syn_info', ...
            'fval_opt', 'ef_opt', 'out_opt', 'lambda_opt', 'grad_opt', 'hess_opt', ...
            'rmse_per_trial', ...
            'sample_list', 'trial_list', 'speed_list', 'leg_list', 'subject_id', 'n');

        % ---- Save prediction comparison figure --------------------------
        fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0 0 1 1]);
        plot_qp_prediction_comparison(subject_id, speed_list, leg_list, ...
            data, rmse_per_trial, Fout(:, :, :, :, 1, 1));
        saveas(fig, [fname_base '.png']);
        close(fig);

        % ---- Save Q / l weight heatmap ----------------------------------
        fig_weights = qp_visualize_weights(Q_opt, l_opt, muscle_names);
        exportgraphics(fig_weights, [fname_base '-weights.png'], 'Resolution', 150);
        close(fig_weights);

        % ---- Save eigendecomposition heatmap ----------------------------
        fig_eig = qp_visualize_eig(Q_opt, muscle_names);
        exportgraphics(fig_eig, [fname_base '-eig.png'], 'Resolution', 150);
        close(fig_eig);

        fprintf('  Saved: %s\n', fname_base);

    end
    end

end
