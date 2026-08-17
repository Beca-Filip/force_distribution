close all; clear all; clc;

% =========================================================================
% QP inverse-optimal-control  --  local fmincon search
%
% Identifies QP cost weights (Q, l) for each subject / speed / leg
% combination by minimising RMSE between QP-predicted and observed forces.
%
% The search variable is theta = [q; l], where
%   q  (n*(n+1)/2 x 1)  lower-triangular entries of the Cholesky factor L
%                        of Q  (Q = L*L', PSD by construction)
%   l  (n x 1)           linear weight vector
%
% Results are saved to main/ (named after this script) as
%   subject-{S}-speed-{speed}-leg-{leg}.mat   (weights + metadata)
%   subject-{S}-speed-{speed}-leg-{leg}.png   (prediction comparison figure)
%
% Each .mat file is self-contained for cross-subject / cross-condition
% validation: load the file, use Q_opt + l_opt in QP_subroutine on any
% other subject or condition.
% =========================================================================

% ---- Subject definitions ------------------------------------------------
subjects(1).id         = 4;
subjects(1).data_dir   = fullfile('..', 'Optimization Model Data', 'Patient4.mat');
subjects(1).speed_list = 1:5;

subjects(2).id         = 5;
subjects(2).data_dir   = fullfile('..', 'Optimization Model Data', 'Patient5.mat');
subjects(2).speed_list = 1:5;

% ---- Output directory ---------------------------------------------------
results_dir = fullfile(fileparts(mfilename('fullpath')), 'main');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% =========================================================================
% Main loop over subjects
% =========================================================================
for si = 1:length(subjects)

    subject_id = subjects(si).id;
    fprintf('=== Subject %d ===\n', subject_id);

    % ---- Load data ------------------------------------------------------
    load(subjects(si).data_dir);

    % Required fields in data for QP:
    %   data.f      [n x 1 x nsamples x ntrials x nspeeds x nlegs]
    %   data.fmin   [n x 1 x nsamples x ntrials x nspeeds x nlegs]
    %   data.fmax   [n x 1 x nsamples x ntrials x nspeeds x nlegs]
    %   data.A      [ne x n x nsamples x ntrials x nspeeds x nlegs]
    %   data.b      [ne x 1 x nsamples x ntrials x nspeeds x nlegs]
    %
    % Fields present in the original DO model that are NOT used here:
    %   data.J_min, data.J_max  (cost normalisation for the DO formulation)
    %   data.vmt, data.M, data.fpassive, data.pcsa, data.f0, data.mass, data.r

    % ---- Bound repair + force-space normalisation -----------------------
    % Shared with the pooled scripts.  It has to be the same code in both, or
    % the per-condition vs pooled comparison compares data preparation as
    % much as it compares objectives.
    data = prepare_qp_data(data);

    % ---- Build QP model -------------------------------------------------
    n = size(data.f, 1);
    [model, vars] = form_casadi_qp_model(n);

    % OSQP tolerances must be TIGHT.  At the defaults the QP solutions are
    % accurate enough for the objective but not for the KKT sensitivity that
    % supplies the gradient, and fmincon then fails its line search and exits
    % with flag -2 ("converged to an infeasible point") after ~20 iterations.
    % Measured on subject 4, speed 1, leg 1: defaults gave exitflag -2 at 20
    % iterations; with these tolerances the same run reached 150 iterations
    % with exitflag 0 and constraint violation 1.6e-3.
    sol_opt = struct;
    sol_opt.osqp.verbose  = 0;
    sol_opt.osqp.eps_abs  = 1e-10;
    sol_opt.osqp.eps_rel  = 1e-10;
    sol_opt.osqp.max_iter = 200000;
    model.solver('osqp', sol_opt);

    % ---- Lists ----------------------------------------------------------
    trial_list  = 1:10;
    sample_list = 1:101;

    % ---- Conditioning ---------------------------------------------------
    % Q = L*L' + (n/cond_max)*I with trace(Q) = n, giving cond(Q) <= cond_max
    % exactly at every feasible iterate.  See theta_to_Ql.m.
    cond_max  = 1e4;
    eps_shift = n / cond_max;

    % ---- Initial theta: Q = I, l = random, feasible for trace(Q) = n ----
    % L0 = sqrt(1-eps_shift)*I  =>  Q0 = (1-eps_shift)*I + eps_shift*I = I,
    % so trace(Q0) = n is satisfied exactly and fmincon starts feasible.
    nq   = n*(n+1)/2;
    mask = tril(true(n));
    [r, c] = find(mask);
    q0   = zeros(nq, 1);
    q0(r == c) = sqrt(1 - eps_shift);
    l0   = randn(n, 1);
    theta0 = [q0; l0];

    % =========================================================================
    % Loop over speeds and legs
    % =========================================================================
    for speed_list = subjects(si).speed_list
    for leg_list = [1, 2]

        fprintf('  Subject %d  Speed %d  Leg %d\n', subject_id, speed_list, leg_list);

        % ---- IO search --------------------------------------------------
        % Timed, because funcCount alone does not tell us how long a fit took
        % and the per-condition runs are the budget estimate for the pooled
        % ones.  summarize_main_results reads t_fit_s.
        t_fit = tic;
        [theta_opt, fval_opt, ef_opt, out_opt, lambda_opt, grad_opt, hess_opt] = ...
            QP_IO_fmincon_search(theta0, data, vars, model, ...
                sample_list, trial_list, speed_list, leg_list, cond_max);
        t_fit_s = toc(t_fit);

        % ---- Reconstruct Q and l ----------------------------------------
        % Must go through theta_to_Ql with the SAME cond_max, otherwise the
        % saved Q_opt would be missing the eps_shift and would not be the
        % matrix that was actually optimised.
        [Q_opt, l_opt, L_opt] = theta_to_Ql(theta_opt, n, cond_max);

        % ---- Conditioning diagnostics -----------------------------------
        eig_Q_opt  = sort(eig(Q_opt));
        cond_Q_opt = eig_Q_opt(end) / eig_Q_opt(1);
        fprintf('    cond(Q) = %.3e   lambda_min = %.3e   trace(Q) = %.4f (n = %d)\n', ...
            cond_Q_opt, eig_Q_opt(1), trace(Q_opt), n);
        if cond_Q_opt > cond_max * (1 + 1e-6)
            warning('main:condExceeded', ...
                'cond(Q) = %.3e exceeds cond_max = %.3e; check constraint feasibility.', ...
                cond_Q_opt, cond_max);
        end

        % ---- Predicted forces (needed for RMSE and figure) --------------
        Fout = QP_subroutine(Q_opt, l_opt, data, vars, model, ...
            sample_list, trial_list, speed_list, leg_list);
        % Fout: [n x 1 x nsamples x ntrials x 1 x 1]

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
            'theta_opt', 'Q_opt', 'L_opt', 'l_opt', ...
            'cond_max', 'eps_shift', 'cond_Q_opt', 'eig_Q_opt', ...
            'fval_opt', 'ef_opt', 'out_opt', 'lambda_opt', 'grad_opt', 'hess_opt', ...
            'rmse_per_trial', 't_fit_s', ...
            'sample_list', 'trial_list', 'speed_list', 'leg_list', 'subject_id', 'n');

        % ---- Save prediction comparison figure --------------------------
        fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0 0 1 1]);
        plot_qp_prediction_comparison(subject_id, speed_list, leg_list, ...
            data, rmse_per_trial, Fout(:, :, :, :, 1, 1));
        saveas(fig, [fname_base '.png']);
        close(fig);

        % ---- Save Q / l weight heatmap ----------------------------------
        fig_weights = qp_visualize_weights(Q_opt, l_opt);
        exportgraphics(fig_weights, [fname_base '-weights.png'], 'Resolution', 150);
        close(fig_weights);

        % ---- Save eigendecomposition heatmap ----------------------------
        fig_eig = qp_visualize_eig(Q_opt);
        exportgraphics(fig_eig, [fname_base '-eig.png'], 'Resolution', 150);
        close(fig_eig);

        fprintf('  Saved: %s\n', fname_base);

    end
    end

end
