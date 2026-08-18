close all; clear all; clc;

% =========================================================================
% QP inverse-optimal-control  --  per-condition fmincon search  (task B5)
%
% Identifies QP cost weights (Q, l) separately for each subject / speed / leg
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
%
%
% SURVIVING A LONG RUN
% --------------------
% Twenty fits at hours apiece, on a machine reached over remote desktop.  Two
% ways that has already been lost, and what is done about each:
%
%   A CasADi error ended the batch at the second condition.  Every QP now
%   goes through QP_SOLVE, which retries a hard QP at a relaxed tolerance;
%   and each condition is wrapped in try/catch, so a condition that cannot be
%   fitted is logged and skipped rather than taking the other nineteen with
%   it.
%
%   The session was lost when the remote desktop rebooted.  QP_LOG writes
%   log.txt, events.jsonl, iterations.csv, status.txt and checkpoint.mat to
%   the results directory, flushing on every write, and RESUME below skips
%   conditions whose .mat is already on disk.  Rerunning after a reboot
%   continues rather than restarts.
%
% Before committing to a full batch, run
%     qp_health_check('../Optimization Model Data/Patient4.mat')
% which answers in about a minute what this script answers in about a day.
% =========================================================================

% ---- Run options --------------------------------------------------------
resume         = true;   % skip conditions whose .mat already exists
seed           = 0;      % l0 is random; without this the run is not reproducible
max_iterations = 1e3;
cond_max       = 1e4;
trial_list     = 1:10;
sample_list    = 1:101;

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

% ---- Logger -------------------------------------------------------------
% Append, because a resumed run must not roll away the log of the run it is
% resuming -- that log is where the failure is recorded.
lg = qp_log(results_dir, 'Append', resume);
cleanup = onCleanup(@() lg.close());

lg.header('B5 -- per-condition QP inverse optimal control', struct( ...
    'resume',         resume, ...
    'seed',           seed, ...
    'max_iterations', max_iterations, ...
    'cond_max',       cond_max, ...
    'n_trials',       numel(trial_list), ...
    'n_samples',      numel(sample_list)));

n_done = 0; n_skipped = 0; n_failed = 0;
t_batch = tic;

% =========================================================================
% Main loop over subjects
% =========================================================================
for si = 1:length(subjects)

    subject_id = subjects(si).id;
    lg.section('Subject %d', subject_id);

    % ---- Load data ------------------------------------------------------
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
    S = load(subjects(si).data_dir);

    % ---- Bound repair + force-space normalisation -----------------------
    % Shared with the pooled scripts.  It has to be the same code in both, or
    % the per-condition vs pooled comparison compares data preparation as
    % much as it compares objectives.
    data = prepare_qp_data(S.data);

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
    %
    % The cost of asking for 1e-10 is that OSQP cannot always deliver it.  It
    % is a first-order method and the Hessian it factorises is W'*Q*W, whose
    % condition number reaches ~1e7 after a hundred fmincon iterations and can
    % reach ~4e10 at the cond(Q) ceiling.  Passing sol_opt to the search below
    % is what makes that survivable: QP_SOLVE relaxes the request for the few
    % QPs that need it, and records having done so, instead of raising.
    % error_on_fail = false is not about ignoring failures -- QP_SOLVE handles
    % them either way.  It changes WHICH layer raises.  Left true, the conic
    % interface raises first and prints every solver input (the full H, A, all
    % bounds: hundreds of lines) to stdout before the message, which in a
    % logged batch means a screenful of numbers per failed QP and a status
    % string that says only "conic process failed".  Set false, Opti raises
    % instead, quietly, and the message carries OSQP's actual return status --
    % which is the one thing worth recording.
    sol_opt = struct;
    sol_opt.error_on_fail = false;
    sol_opt.osqp.verbose  = 0;
    sol_opt.osqp.eps_abs  = 1e-10;
    sol_opt.osqp.eps_rel  = 1e-10;
    sol_opt.osqp.max_iter = 200000;
    model.solver('osqp', sol_opt);

    % ---- Conditioning ---------------------------------------------------
    % Q = L*L' + (n/cond_max)*I with trace(Q) = n, giving cond(Q) <= cond_max
    % exactly at every feasible iterate.  See theta_to_Ql.m.
    eps_shift = n / cond_max;

    % ---- Initial theta: Q = I, l = random, feasible for trace(Q) = n ----
    % L0 = sqrt(1-eps_shift)*I  =>  Q0 = (1-eps_shift)*I + eps_shift*I = I,
    % so trace(Q0) = n is satisfied exactly and fmincon starts feasible.
    %
    % The seed matters.  l0 was previously an unseeded randn, which made a
    % rerun a different experiment from the run it was meant to reproduce.
    nq   = n*(n+1)/2;
    mask = tril(true(n));
    [r, c] = find(mask);
    q0   = zeros(nq, 1);
    q0(r == c) = sqrt(1 - eps_shift);
    rng(seed, 'twister');
    l0   = randn(n, 1);
    theta0 = [q0; l0];

    lg.printf('n = %d muscles, eps_shift = %.4e, cond(W) = %.3e', ...
        n, eps_shift, cond(data.F_invcov));

    % =====================================================================
    % Loop over speeds and legs
    % =====================================================================
    for speed_list = subjects(si).speed_list
    for leg_list = [1, 2]

        tag = sprintf('subject-%d-speed-%d-leg-%d', subject_id, speed_list, leg_list);
        fname_base = fullfile(results_dir, tag);

        if resume && exist([fname_base '.mat'], 'file')
            lg.printf('SKIP %s (already on disk)', tag);
            n_skipped = n_skipped + 1;
            continue;
        end

        lg.section('Subject %d  Speed %d  Leg %d', subject_id, speed_list, leg_list);
        lg.status('%s: starting', tag);
        lg.event('condition_start', struct('subject', subject_id, ...
            'speed', speed_list, 'leg', leg_list));

        try
            % ---- IO search ----------------------------------------------
            % Timed, because funcCount alone does not tell us how long a fit
            % took and the per-condition runs are the budget estimate for the
            % pooled ones.  summarize_main_results reads t_fit_s.
            fmc_overrides = struct( ...
                'PlotFcn',       {{}}, ...   % batch run: no figure windows
                'Display',       'iter', ...
                'MaxIterations', max_iterations, ...
                'OutputFcn',     @(th, ov, st) lg.fmincon_hook(th, ov, st, tag));

            t_fit = tic;
            [theta_opt, fval_opt, ef_opt, out_opt, lambda_opt, grad_opt, hess_opt] = ...
                QP_IO_fmincon_search(theta0, data, vars, model, ...
                    sample_list, trial_list, speed_list, leg_list, cond_max, ...
                    fmc_overrides, sol_opt);
            t_fit_s = toc(t_fit);

            % ---- Reconstruct Q and l ------------------------------------
            % Must go through theta_to_Ql with the SAME cond_max, otherwise
            % the saved Q_opt would be missing the eps_shift and would not be
            % the matrix that was actually optimised.
            [Q_opt, l_opt, L_opt] = theta_to_Ql(theta_opt, n, cond_max);

            % ---- Conditioning diagnostics -------------------------------
            eig_Q_opt  = sort(eig(Q_opt));
            cond_Q_opt = eig_Q_opt(end) / eig_Q_opt(1);
            H_opt      = data.F_invcov' * Q_opt * data.F_invcov;
            cond_H_opt = cond((H_opt + H_opt') / 2);

            lg.printf(['done in %s: exitflag %d, fval %.4f, %d iterations, ' ...
                'firstorderopt %.3e'], hms(t_fit_s), ef_opt, fval_opt, ...
                out_opt.iterations, out_opt.firstorderopt);
            lg.printf('cond(Q) = %.3e   lambda_min = %.3e   trace(Q) = %.4f (n = %d)', ...
                cond_Q_opt, eig_Q_opt(1), trace(Q_opt), n);
            lg.printf('cond(W''*Q*W) = %.3e   <- what OSQP had to factorise', cond_H_opt);

            if cond_Q_opt > cond_max * (1 + 1e-6)
                lg.warn('main:condExceeded', ...
                    'cond(Q) = %.3e exceeds cond_max = %.3e; check constraint feasibility.', ...
                    cond_Q_opt, cond_max);
            end
            if ef_opt == 2
                lg.warn('main:stalledOnStep', ...
                    ['exit flag 2: the step fell below tolerance while firstorderopt ' ...
                     'was still %.3g.  Expected on this problem -- the pooled objective ' ...
                     'is only piecewise smooth -- but it means the fit stopped, not converged.'], ...
                    out_opt.firstorderopt);
            end

            % ---- Predicted forces (needed for RMSE and figure) ----------
            Fout = QP_subroutine(Q_opt, l_opt, data, vars, model, ...
                sample_list, trial_list, speed_list, leg_list, [], sol_opt);
            % Fout: [n x 1 x nsamples x ntrials x 1 x 1]

            % ---- Per-trial RMSE -----------------------------------------
            ntrials = length(trial_list);
            rmse_per_trial = zeros(ntrials, 1);
            for ti = 1:ntrials
                f_ref_t  = data.f(:, :, sample_list, trial_list(ti), speed_list, leg_list);
                f_pred_t = Fout(:, :, :, ti, 1, 1);
                rmse_per_trial(ti) = rmse(f_ref_t(:), f_pred_t(:));
            end

            % ---- Save .mat ----------------------------------------------
            save([fname_base '.mat'], ...
                'theta_opt', 'Q_opt', 'L_opt', 'l_opt', ...
                'cond_max', 'eps_shift', 'cond_Q_opt', 'eig_Q_opt', 'cond_H_opt', ...
                'fval_opt', 'ef_opt', 'out_opt', 'lambda_opt', 'grad_opt', 'hess_opt', ...
                'rmse_per_trial', 't_fit_s', 'seed', ...
                'sample_list', 'trial_list', 'speed_list', 'leg_list', 'subject_id', 'n');

            % ---- Save prediction comparison figure ----------------------
            fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0 0 1 1]);
            plot_qp_prediction_comparison(subject_id, speed_list, leg_list, ...
                data, rmse_per_trial, Fout(:, :, :, :, 1, 1));
            saveas(fig, [fname_base '.png']);
            close(fig);

            % ---- Save Q / l weight heatmap ------------------------------
            fig_weights = qp_visualize_weights(Q_opt, l_opt);
            exportgraphics(fig_weights, [fname_base '-weights.png'], 'Resolution', 150);
            close(fig_weights);

            % ---- Save eigendecomposition heatmap ------------------------
            fig_eig = qp_visualize_eig(Q_opt);
            exportgraphics(fig_eig, [fname_base '-eig.png'], 'Resolution', 150);
            close(fig_eig);

            lg.printf('saved %s', fname_base);
            lg.event('condition_done', struct( ...
                'subject', subject_id, 'speed', speed_list, 'leg', leg_list, ...
                'fval', fval_opt, 'exitflag', ef_opt, ...
                'iterations', out_opt.iterations, ...
                'firstorderopt', out_opt.firstorderopt, ...
                'cond_Q', cond_Q_opt, 'cond_H', cond_H_opt, ...
                't_fit_s', t_fit_s));
            n_done = n_done + 1;

        catch ME
            % One condition failing is not a reason to lose the other
            % nineteen.  Everything needed to reproduce it in isolation is in
            % the log; the batch moves on.
            lg.exception(ME, struct('subject', subject_id, ...
                'speed', speed_list, 'leg', leg_list, 'tag', tag));
            lg.printf('CONTINUING with the next condition.');
            n_failed = n_failed + 1;
            close all force;
        end

    end
    end

end

lg.printf('batch finished in %s: %d fitted, %d skipped, %d failed', ...
    hms(toc(t_batch)), n_done, n_skipped, n_failed);
if n_failed > 0
    lg.printf('rerun this script to retry the failed conditions (resume skips the rest).');
end

% =========================================================================
function s = hms(t)
s = sprintf('%02d:%02d:%02d', floor(t/3600), mod(floor(t/60), 60), mod(floor(t), 60));
end
