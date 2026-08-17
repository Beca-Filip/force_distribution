function summary = run_pooled_fit(cfg)
%RUN_POOLED_FIT  Fit ONE (Q, l) jointly across all conditions of one subject.
%
%   summary = RUN_POOLED_FIT(cfg)
%
%   This is the engine behind main_pooled_s1.m and main_pooled_s2.m.  Where
%   main.m runs one fmincon per (speed, leg) and produces 10 different cost
%   functions per subject, this runs ONE fmincon over every (sample, trial,
%   speed, leg) at once and produces a single cost function that has to
%   explain all of them.  That is the whole question: can one QP do it?
%
%   The multi-condition objective and gradient are already pooled inside
%   QP_IO_INNER_LOOP, so the only change here is passing vector speed_list and
%   leg_list instead of scalars.
%
%   COST.  One objective evaluation solves
%   numel(sample_list)*numel(trial_list)*numel(speed_list)*numel(leg_list)
%   QPs -- for the full lists that is 10,100 against 1,010 for a per-condition
%   fit, so a pooled fit costs roughly 10x one entry of the B5 batch at equal
%   iteration counts.  Budget accordingly, and note that cfg.max_iterations is
%   the knob that actually bounds the runtime.
%
%   CHECKPOINTING.  Because the run is long, an OutputFcn writes
%   <results_dir>/checkpoint.mat after every iteration (iteration number,
%   theta, fval, constraint violation, elapsed time).  An interrupted run can
%   be restarted from it by passing cfg.theta0.
%
%   OUTPUT LAYOUT.  Everything lands in cfg.results_dir:
%     fit.mat                          the pooled fmincon record
%     fit-weights.png, fit-eig.png     Q and its eigendecomposition (one Q,
%                                      so these are folder-level, not per
%                                      condition)
%     subject-S-speed-P-leg-L.mat      the pooled (Q,l) SCORED on one
%                                      condition, one file per condition
%     subject-S-speed-P-leg-L.png      prediction comparison for it
%     checkpoint.mat, log.txt
%
%   The per-condition .mat files deliberately use the same variable names as
%   main.m's, so SUMMARIZE_MAIN_RESULTS runs on this folder unchanged and the
%   per-condition vs pooled comparison is a table join.  Two differences are
%   flagged inside those files:
%     is_pooled_eval = true   they are evaluations, not independent fits
%     fval_opt                the objective RESTRICTED to this condition;
%                             fval_pooled carries the joint objective
%     t_fit_s                 the time to SCORE this condition;
%                             t_pooled_fit_s carries the fit time
%   Consequently SUMMARIZE_MAIN_RESULTS's "saved (Q,l) reproduces fval" check
%   is vacuous on this folder -- fval_opt is computed from the same solve.
%   The equivalent check with teeth is run here instead: the quadratic mean of
%   the per-condition RMSEs must reproduce the pooled objective fmincon
%   returned, and it is a genuine re-solve.
%
%   CFG FIELDS
%     subject_id      integer, used in file names and titles
%     data_path       path to the Patient*.mat
%     results_dir     output directory (created if absent)
%   optional:
%     speed_list      default 1:size(data.f,5)
%     leg_list        default 1:size(data.f,6)
%     sample_list     default 1:101
%     trial_list      default 1:10
%     cond_max        default 1e4
%     max_iterations  default 300  (fmincon's own default of 1000 is not
%                     affordable here; see COST above)
%     seed            default 0, for the random l0
%     theta0          default [], meaning Q0 = I and l0 = randn -- pass a
%                     checkpoint's theta to resume
%     make_figures    default true
%
%   See also MAIN, QP_IO_FMINCON_SEARCH, QP_IO_INNER_LOOP, PREPARE_QP_DATA,
%   SUMMARIZE_MAIN_RESULTS.

% =========================================================================
% Config
% =========================================================================
required = {'subject_id', 'data_path', 'results_dir'};
for k = 1:numel(required)
    if ~isfield(cfg, required{k})
        error('run_pooled_fit:missingConfig', 'cfg.%s is required.', required{k});
    end
end

cfg = default_field(cfg, 'sample_list',    1:101);
cfg = default_field(cfg, 'trial_list',     1:10);
cfg = default_field(cfg, 'cond_max',       1e4);
cfg = default_field(cfg, 'max_iterations', 300);
cfg = default_field(cfg, 'seed',           0);
cfg = default_field(cfg, 'theta0',         []);
cfg = default_field(cfg, 'make_figures',   true);

assert_environment();

results_dir = cfg.results_dir;
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

log_path = fullfile(results_dir, 'log.txt');
diary(log_path);
diary on;
cleanup = onCleanup(@() diary('off'));

fprintf('=====================================================================\n');
fprintf(' POOLED FIT -- subject %d\n', cfg.subject_id);
fprintf(' started %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(' output  %s\n', results_dir);
fprintf('=====================================================================\n');

% =========================================================================
% Data
% =========================================================================
S    = load(cfg.data_path);
data = prepare_qp_data(S.data);

n = size(data.f, 1);

cfg = default_field(cfg, 'speed_list', 1:size(data.f, 5));
cfg = default_field(cfg, 'leg_list',   1:size(data.f, 6));

sample_list = cfg.sample_list;
trial_list  = cfg.trial_list;
speed_list  = cfg.speed_list;
leg_list    = cfg.leg_list;

check_range('speed_list', speed_list, size(data.f, 5));
check_range('leg_list',   leg_list,   size(data.f, 6));
check_range('sample_list', sample_list, size(data.f, 3));
check_range('trial_list',  trial_list,  size(data.f, 4));

n_qp_per_eval = numel(sample_list) * numel(trial_list) * ...
                numel(speed_list)  * numel(leg_list);

fprintf('\n  muscles      : %d\n', n);
fprintf('  speeds       : %s\n', mat2str(speed_list));
fprintf('  legs         : %s\n', mat2str(leg_list));
fprintf('  samples      : %d\n', numel(sample_list));
fprintf('  trials       : %d\n', numel(trial_list));
fprintf('  QPs per objective evaluation : %d\n', n_qp_per_eval);
fprintf('  cond_max     : %g\n', cfg.cond_max);
fprintf('  max iterations : %d\n\n', cfg.max_iterations);

% =========================================================================
% Model
% =========================================================================
[model, vars] = form_casadi_qp_model(n);

% Tight OSQP, for the reason documented in main.m: at the defaults the QP
% solutions are accurate enough for the objective but not for the KKT
% sensitivity that supplies the gradient, and fmincon's line search then
% stalls.  This matters more here, not less -- the pooled gradient is a sum
% over ten times as many solves, so ten times as much solver noise.
sol_opt = struct;
sol_opt.osqp.verbose  = 0;
sol_opt.osqp.eps_abs  = 1e-10;
sol_opt.osqp.eps_rel  = 1e-10;
sol_opt.osqp.max_iter = 200000;
model.solver('osqp', sol_opt);

% =========================================================================
% Initial theta
% =========================================================================
% Q0 = I, exactly as in main.m, and NOT warm-started from the per-condition
% fits.  Same start means any difference between the pooled and per-condition
% results is attributable to pooling rather than to initialisation.
cond_max  = cfg.cond_max;
eps_shift = n / cond_max;

if isempty(cfg.theta0)
    rng(cfg.seed, 'twister');   % l0 is random; fix it so a long run is repeatable
    nq   = n*(n+1)/2;
    mask = tril(true(n));
    [r, c] = find(mask);
    q0   = zeros(nq, 1);
    q0(r == c) = sqrt(1 - eps_shift);   % L0 = sqrt(1-eps)*I => Q0 = I, trace = n
    l0   = randn(n, 1);
    theta0 = [q0; l0];
    fprintf('  theta0: Q0 = I, l0 = randn (seed %d)\n', cfg.seed);
else
    theta0 = reshape(cfg.theta0, [], 1);
    fprintf('  theta0: supplied (resuming)\n');
end

% =========================================================================
% Fit
% =========================================================================
ckpt_path = fullfile(results_dir, 'checkpoint.mat');
t_fit     = tic;

fmc_overrides = struct( ...
    'PlotFcn',       {{}}, ...            % batch run: no figure windows
    'Display',       'iter', ...
    'MaxIterations', cfg.max_iterations, ...
    'OutputFcn',     @(th, ov, st) checkpoint_fcn(th, ov, st, ckpt_path, t_fit));

fprintf('\n--- fmincon ---------------------------------------------------------\n');
[theta_opt, fval_opt, ef_opt, out_opt, lambda_opt, grad_opt, hess_opt] = ...
    QP_IO_fmincon_search(theta0, data, vars, model, ...
        sample_list, trial_list, speed_list, leg_list, cond_max, fmc_overrides);

t_pooled_fit_s = toc(t_fit);

[Q_opt, l_opt, L_opt] = theta_to_Ql(theta_opt, n, cond_max);

eig_Q_opt  = sort(eig((Q_opt + Q_opt') / 2));
cond_Q_opt = eig_Q_opt(end) / eig_Q_opt(1);

fprintf('\n--- fit --------------------------------------------------------------\n');
fprintf('  exitflag    : %d\n', ef_opt);
fprintf('  iterations  : %d   funcCount %d\n', out_opt.iterations, out_opt.funcCount);
fprintf('  pooled E    : %.6f\n', fval_opt);
fprintf('  time        : %.1f s (%.2f h)\n', t_pooled_fit_s, t_pooled_fit_s/3600);
fprintf('  cond(Q)     : %.3e   lambda_min %.3e   trace(Q) %.4f (n = %d)\n', ...
    cond_Q_opt, eig_Q_opt(1), trace(Q_opt), n);

if cond_Q_opt > cond_max * (1 + 1e-6)
    warning('run_pooled_fit:condExceeded', ...
        'cond(Q) = %.3e exceeds cond_max = %.3e; check constraint feasibility.', ...
        cond_Q_opt, cond_max);
end

subject_id = cfg.subject_id;

save(fullfile(results_dir, 'fit.mat'), ...
    'theta_opt', 'Q_opt', 'L_opt', 'l_opt', ...
    'cond_max', 'eps_shift', 'cond_Q_opt', 'eig_Q_opt', ...
    'fval_opt', 'ef_opt', 'out_opt', 'lambda_opt', 'grad_opt', 'hess_opt', ...
    't_pooled_fit_s', 'theta0', ...
    'sample_list', 'trial_list', 'speed_list', 'leg_list', ...
    'subject_id', 'n');

% =========================================================================
% Score the pooled (Q, l) on each condition separately
% =========================================================================
% Per-condition RMSE is what makes the pooled fit interpretable: a single
% number over ten conditions hides whether the QP explains them evenly or
% explains two of them and fails the rest.
fprintf('\n--- per-condition evaluation ----------------------------------------\n');

n_cond    = numel(speed_list) * numel(leg_list);
cond_rows = cell(n_cond, 1);
ci        = 0;

for speed = speed_list
for leg   = leg_list

    ci = ci + 1;
    t_eval = tic;

    Fout = QP_subroutine(Q_opt, l_opt, data, vars, model, ...
        sample_list, trial_list, speed, leg);      % [n x 1 x nsamp x ntrial x 1 x 1]

    ntrials        = numel(trial_list);
    rmse_per_trial = zeros(ntrials, 1);
    for ti = 1:ntrials
        f_ref_t  = data.f(:, :, sample_list, trial_list(ti), speed, leg);
        f_pred_t = Fout(:, :, :, ti, 1, 1);
        rmse_per_trial(ti) = rmse(f_ref_t(:), f_pred_t(:));
    end

    f_ref_c  = data.f(:, :, sample_list, trial_list, speed, leg);
    rmse_cond = rmse(f_ref_c(:), reshape(Fout, [], 1));

    t_fit_s = toc(t_eval);

    fprintf('  speed %d leg %d : E_cond = %8.4f   (trials %.3f to %.3f)   %.1f s\n', ...
        speed, leg, rmse_cond, min(rmse_per_trial), max(rmse_per_trial), t_fit_s);

    % ---- save, using main.m's variable names -----------------------------
    fname_base = fullfile(results_dir, ...
        sprintf('subject-%d-speed-%d-leg-%d', cfg.subject_id, speed, leg));

    % Assembled as a struct and written with -struct, so that fval_opt,
    % speed_list and leg_list carry THIS condition's values while the
    % workspace variables of the same name keep the pooled ones.
    % grad_opt and hess_opt are not copied: they belong to the pooled fit and
    % are 665 and 665x665 respectively.
    rec = struct();
    rec.theta_opt         = theta_opt;
    rec.Q_opt             = Q_opt;
    rec.L_opt             = L_opt;
    rec.l_opt             = l_opt;
    rec.cond_max          = cond_max;
    rec.eps_shift         = eps_shift;
    rec.cond_Q_opt        = cond_Q_opt;
    rec.eig_Q_opt         = eig_Q_opt;
    rec.fval_opt          = rmse_cond;      % this condition
    rec.fval_pooled       = fval_opt;       % the joint objective
    rec.ef_opt            = ef_opt;
    rec.out_opt           = out_opt;
    rec.lambda_opt        = lambda_opt;
    rec.rmse_per_trial    = rmse_per_trial;
    rec.t_fit_s           = t_fit_s;        % time to SCORE this condition
    rec.t_pooled_fit_s    = t_pooled_fit_s; % time to FIT
    rec.sample_list       = sample_list;
    rec.trial_list        = trial_list;
    rec.speed_list        = speed;
    rec.leg_list          = leg;
    rec.subject_id        = cfg.subject_id;
    rec.n                 = n;
    rec.is_pooled_eval    = true;
    rec.pooled_speed_list = speed_list;
    rec.pooled_leg_list   = leg_list;

    save([fname_base '.mat'], '-struct', 'rec');

    if cfg.make_figures
        fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0 0 1 1]);
        plot_qp_prediction_comparison(cfg.subject_id, speed, leg, ...
            data, rmse_per_trial, Fout(:, :, :, :, 1, 1));
        saveas(fig, [fname_base '.png']);
        close(fig);
    end

    cond_rows{ci} = struct('speed', speed, 'leg', leg, ...
        'rmse_cond', rmse_cond, ...
        'rmse_trial_min', min(rmse_per_trial), ...
        'rmse_trial_max', max(rmse_per_trial), ...
        'rmse_trial_std', std(rmse_per_trial), ...
        't_eval_s', t_fit_s);
end
end

cond_table = struct2table([cond_rows{:}]);

% =========================================================================
% Consistency check with teeth
% =========================================================================
% Every condition contributes the same number of residuals, so the RMSE over
% all of them is the quadratic mean of the per-condition RMSEs.  These
% evaluations are independent re-solves, so if this does not reproduce the
% objective fmincon reported, the QP minimiser is not reproducible at this Q
% and the fit means nothing -- the same failure four of the legacy
% per-condition fits showed.
rmse_recomputed = sqrt(mean(cond_table.rmse_cond .^ 2));
pooled_resid    = abs(rmse_recomputed - fval_opt);
pooled_ok       = pooled_resid < 1e-6 * max(1, fval_opt);

fprintf('\n--- consistency -----------------------------------------------------\n');
fprintf('  fmincon pooled E            : %.6f\n', fval_opt);
fprintf('  re-solved, quadratic mean   : %.6f\n', rmse_recomputed);
if pooled_ok
    fprintf('  [PASS] the pooled objective is reproducible (|diff| %.2e)\n', pooled_resid);
else
    fprintf('  [FAIL] |diff| = %.3e -- re-solving the QP with the saved (Q,l) does\n', pooled_resid);
    fprintf('         not reproduce the objective. The minimiser is not unique at\n');
    fprintf('         this Q, or OSQP is not converging. Do not use this fit.\n');
end

% =========================================================================
% Weight figures (one Q, so folder-level rather than per condition)
% =========================================================================
if cfg.make_figures
    fig_weights = qp_visualize_weights(Q_opt, l_opt);
    exportgraphics(fig_weights, fullfile(results_dir, 'fit-weights.png'), 'Resolution', 150);
    close(fig_weights);

    fig_eig = qp_visualize_eig(Q_opt);
    exportgraphics(fig_eig, fullfile(results_dir, 'fit-eig.png'), 'Resolution', 150);
    close(fig_eig);
end

% =========================================================================
% Summary
% =========================================================================
summary = struct();
summary.subject_id      = cfg.subject_id;
summary.n               = n;
summary.results_dir     = results_dir;
summary.Q_opt           = Q_opt;
summary.l_opt           = l_opt;
summary.theta_opt       = theta_opt;
summary.fval_pooled     = fval_opt;
summary.rmse_recomputed = rmse_recomputed;
summary.pooled_ok       = pooled_ok;
summary.exitflag        = ef_opt;
summary.iterations      = out_opt.iterations;
summary.funcCount       = out_opt.funcCount;
summary.firstorderopt   = out_opt.firstorderopt;
summary.constrviolation = out_opt.constrviolation;
summary.cond_Q          = cond_Q_opt;
summary.lambda_min      = eig_Q_opt(1);
summary.trace_Q         = trace(Q_opt);
summary.t_pooled_fit_s  = t_pooled_fit_s;
summary.cond_table      = cond_table;

save(fullfile(results_dir, 'summary.mat'), '-struct', 'summary');

fprintf('\n--- per-condition summary -------------------------------------------\n');
disp(cond_table);
fprintf('  pooled E %.4f over %d conditions;  worst condition %.4f, best %.4f\n', ...
    fval_opt, n_cond, max(cond_table.rmse_cond), min(cond_table.rmse_cond));
fprintf('  finished %s   total fit time %.2f h\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), t_pooled_fit_s/3600);
fprintf('=====================================================================\n');

diary off;

end

% =========================================================================
% Helpers
% =========================================================================

function assert_environment()
%ASSERT_ENVIRONMENT  Fail loudly now rather than wrongly in six hours.
%
%   Checked here and not in main.m because a pooled fit is a long unattended
%   run, usually on a different machine, and two of these failures are silent.

this_dir = fileparts(mfilename('fullpath'));   % workspace_qp
repo_dir = fileparts(this_dir);

addpath(this_dir);
addpath(fullfile(repo_dir, 'workspace_do'));           % eq/ineq_constraint_function
addpath(fullfile(repo_dir, 'utils', 'error_utils'));   % scalar rmse

% CasADi is not in the repo, so it can only be found or reported missing.
if exist('casadi.Opti', 'class') ~= 8
    casadi_path = 'C:/Users/filip/Documents/GitHub/casadi-3.7.0-windows64-matlab2018b';
    if exist(casadi_path, 'dir')
        addpath(casadi_path);
    end
end
if exist('casadi.Opti', 'class') ~= 8
    error('run_pooled_fit:noCasadi', ...
        'CasADi is not on the path. addpath your casadi-3.7.0 directory first.');
end

% The project rmse flattens with a(:) and returns a SCALAR; MATLAB's R2022b+
% built-in reduces along dim 1 and returns an ARRAY.  With the built-in bound
% the objective is silently non-scalar.  See run_qp_tests.m.
if ~isequal(size(rmse(ones(2,2,2), zeros(2,2,2))), [1 1])
    error('run_pooled_fit:wrongRmse', ...
        ['MATLAB''s built-in rmse is shadowing the project version in %s. ' ...
         'The objective would be an array instead of a scalar.'], ...
        fullfile(repo_dir, 'utils', 'error_utils'));
end

if isempty(which('eq_constraint_function'))
    error('run_pooled_fit:noConstraints', ...
        'eq_constraint_function is not on the path (expected in %s).', ...
        fullfile(repo_dir, 'workspace_do'));
end
end

function cfg = default_field(cfg, name, value)
if ~isfield(cfg, name) || isempty(cfg.(name))
    cfg.(name) = value;
end
end

function check_range(name, list, dim_size)
if any(list < 1) || any(list > dim_size)
    error('run_pooled_fit:badList', ...
        '%s = %s is outside the data range 1..%d.', name, mat2str(list), dim_size);
end
end

function stop = checkpoint_fcn(theta, optimValues, state, ckpt_path, t_fit)
%CHECKPOINT_FCN  Write the current iterate after every fmincon iteration.
%   A pooled fit runs for hours; without this an interruption loses all of it.
%   Resume by passing the saved theta as cfg.theta0.
stop = false;
if ~any(strcmp(state, {'iter', 'done'}))
    return
end

ck = struct();
ck.iteration       = optimValues.iteration;
ck.theta           = theta;
ck.fval            = optimValues.fval;
ck.firstorderopt   = get_or_nan(optimValues, 'firstorderopt');
ck.constrviolation = get_or_nan(optimValues, 'constrviolation');
ck.elapsed_s       = toc(t_fit);
ck.saved_at        = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% A failed checkpoint write must never abort a multi-hour fit.
try
    save(ckpt_path, '-struct', 'ck');
catch err
    warning('run_pooled_fit:checkpointFailed', ...
        'Could not write %s: %s', ckpt_path, err.message);
end
end

function v = get_or_nan(s, name)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = NaN;
end
end
