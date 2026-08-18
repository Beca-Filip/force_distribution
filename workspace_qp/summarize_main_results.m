function T = summarize_main_results(results_dir, varargin)
%SUMMARIZE_MAIN_RESULTS  Post-process a main/ IOC run and audit it.
%
%   T = SUMMARIZE_MAIN_RESULTS()                 uses ./main
%   T = SUMMARIZE_MAIN_RESULTS(dir)
%   T = SUMMARIZE_MAIN_RESULTS(dir, 'Csv', true, 'Plot', true)
%
%   Reads every subject-*-speed-*-leg-*.mat in results_dir and answers three
%   questions, in this order:
%
%     1. PROVENANCE -- which version of main.m produced these files?  The
%        conditioning fix (B) added cond_max / eps_shift / cond_Q_opt /
%        eig_Q_opt to the save list and a nonlinear equality to fmincon.  A
%        file missing those variables, or carrying an empty lambda.eqnonlin,
%        was produced by code that predates the fix.  This check comes first
%        because every number below it is meaningless if it fails.
%
%     2. CONVERGENCE -- exit flags, iteration and function-evaluation counts,
%        first-order optimality, constraint violation, wall-clock time.
%
%     3. SOLUTION QUALITY -- objective, per-trial RMSE spread, and the
%        conditioning of the fitted Q (trace, lambda_min, cond).
%
%   Each check prints PASS / WARN / FAIL.  A FAIL means the fits should not be
%   used; a WARN means look before you trust.
%
%   Wall-clock time is only reported if main.m saved it (variable t_fit_s).
%   Runs that did not save it show NaN, and funcCount is the cost proxy: one
%   function evaluation is one full pass over every sample/trial, each of
%   which solves a QP, so it is close to proportional to runtime.
%
%   Returns a table with one row per fit, sorted by subject, speed, leg.

if nargin < 1 || isempty(results_dir)
    results_dir = fullfile(fileparts(mfilename('fullpath')), 'main');
end

p = inputParser;
p.addParameter('Csv',  false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Plot', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
want_csv  = logical(p.Results.Csv);
want_plot = logical(p.Results.Plot);

files = dir(fullfile(results_dir, 'subject-*-speed-*-leg-*.mat'));
if isempty(files)
    error('summarize_main_results:noFiles', ...
        'No subject-*-speed-*-leg-*.mat found in %s.', results_dir);
end

% =====================================================================
% Collect
% =====================================================================
nf = numel(files);
rows = cell(nf, 1);

for k = 1:nf
    f = fullfile(results_dir, files(k).name);
    S = load(f);

    r = struct();
    r.file    = string(files(k).name);
    r.subject = get_or(S, 'subject_id', NaN);
    r.speed   = get_or(S, 'speed_list', NaN);
    r.leg     = get_or(S, 'leg_list',   NaN);
    r.n       = get_or(S, 'n',          NaN);
    r.nsample = numel(get_or(S, 'sample_list', []));
    r.ntrial  = numel(get_or(S, 'trial_list',  []));

    % ---- provenance -------------------------------------------------
    % Two independent signals: what main.m chose to save, and whether
    % fmincon was actually handed a nonlinear equality.  They can disagree
    % (cond_max = Inf saves the variable but passes nonlcon = []), which is
    % why both are recorded.
    r.has_cond_vars = all(isfield(S, {'cond_max', 'eps_shift', 'cond_Q_opt', 'eig_Q_opt'}));
    r.cond_max      = get_or(S, 'cond_max',  NaN);
    r.eps_shift     = get_or(S, 'eps_shift', NaN);

    lam = get_or(S, 'lambda_opt', struct());
    if isstruct(lam) && isfield(lam, 'eqnonlin')
        r.n_eqnonlin = numel(lam.eqnonlin);
    else
        r.n_eqnonlin = 0;
    end
    r.constrained = r.n_eqnonlin > 0;

    % Written by run_pooled_fit; absent in main/ files, which are independent
    % per-condition fits.
    r.is_pooled_eval = logical(get_or(S, 'is_pooled_eval', false));

    % ---- convergence ------------------------------------------------
    r.exitflag = get_or(S, 'ef_opt', NaN);
    r.exit_msg = string(exitflag_label(r.exitflag));

    out = get_or(S, 'out_opt', struct());
    r.iterations    = get_or(out, 'iterations',      NaN);
    r.funcCount     = get_or(out, 'funcCount',       NaN);
    r.constrviol    = get_or(out, 'constrviolation', NaN);
    r.firstorderopt = get_or(out, 'firstorderopt',   NaN);
    r.stepsize      = get_or(out, 'stepsize',        NaN);
    r.algorithm     = string(get_or(out, 'algorithm', "unknown"));

    % main.m only stores this if it was written to (see the tic/toc there).
    r.time_s = get_or(S, 't_fit_s', NaN);

    % ---- solution ---------------------------------------------------
    r.fval = get_or(S, 'fval_opt', NaN);

    rp = get_or(S, 'rmse_per_trial', []);
    if isempty(rp)
        [r.rmse_mean, r.rmse_std, r.rmse_min, r.rmse_max, r.rmse_pooled] = deal(NaN);
    else
        rp = rp(:);
        r.rmse_mean   = mean(rp);
        r.rmse_std    = std(rp);
        r.rmse_min    = min(rp);
        r.rmse_max    = max(rp);
        % Every trial contributes the same number of residuals, so the RMSE
        % over all of them is the quadratic mean of the per-trial RMSEs.
        % This must reproduce fval: if it does not, the saved (Q,l) is not
        % the point fmincon reported.
        r.rmse_pooled = sqrt(mean(rp.^2));
    end
    r.fval_resid = abs(r.fval - r.rmse_pooled);

    % ---- conditioning -----------------------------------------------
    % Recomputed from the saved Q rather than read from cond_Q_opt, so this
    % works on files that predate those variables and double-checks the ones
    % that have them.
    Q = get_or(S, 'Q_opt', []);
    if isempty(Q)
        [r.trace_Q, r.lambda_min, r.lambda_max, r.cond_Q] = deal(NaN);
    else
        e = sort(eig((Q + Q') / 2));
        r.trace_Q    = trace(Q);
        r.lambda_min = e(1);
        r.lambda_max = e(end);
        % Guard the ratio: a non-positive lambda_min makes cond meaningless,
        % and a negative one silently produces a negative "condition number".
        if e(1) > 0
            r.cond_Q = e(end) / e(1);
        else
            r.cond_Q = Inf;
        end
    end
    r.trace_resid = abs(r.trace_Q - r.n);

    l = get_or(S, 'l_opt', []);
    r.norm_l = norm(l(:));

    rows{k} = r;
end

T = struct2table([rows{:}]);
T = sortrows(T, {'subject', 'speed', 'leg'});

% =====================================================================
% Report
% =====================================================================
fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' main/ run summary -- %s\n', results_dir);
fprintf(' %d fits, %d subject(s), %s\n', height(T), ...
    numel(unique(T.subject)), datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('=====================================================================\n');

n_fail = 0;
n_warn = 0;

% ---------------------------------------------------------------------
fprintf('\n--- 1. PROVENANCE ---------------------------------------------------\n');

if all(T.has_cond_vars)
    verdict('PASS', 'every file records cond_max / eps_shift / cond_Q_opt / eig_Q_opt');
elseif ~any(T.has_cond_vars)
    n_fail = n_fail + 1;
    verdict('FAIL', ['no file records cond_max -- these were produced by a main.m ' ...
        'that predates the conditioning fix (B).']);
else
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf('%d of %d files record cond_max -- the run mixes code versions.', ...
        sum(T.has_cond_vars), height(T)));
end

if all(T.constrained)
    verdict('PASS', 'fmincon carried the trace(Q) = n equality in every fit');
elseif ~any(T.constrained)
    n_fail = n_fail + 1;
    verdict('FAIL', ['lambda.eqnonlin is empty in every fit -- fmincon ran with ' ...
        'nonlcon = [], so the conditioning constraint was never imposed.']);
else
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf('%d of %d fits were constrained.', sum(T.constrained), height(T)));
end

% A pooled folder holds one fit scored on each condition, not independent
% fits.  Say so, because two checks below read differently there.
if all(T.is_pooled_eval)
    verdict('INFO', ['pooled-evaluation folder: these are one (Q,l) scored on each ' ...
        'condition, so the convergence block describes a single fit and the ' ...
        'fval-reproduction check is vacuous (fval_opt comes from the same solve).']);
elseif any(T.is_pooled_eval)
    n_warn = n_warn + 1;
    verdict('WARN', 'the folder mixes pooled evaluations with independent fits.');
end

algs = unique(T.algorithm);
if numel(algs) == 1
    verdict('INFO', sprintf('algorithm: %s', algs(1)));
else
    n_warn = n_warn + 1;
    verdict('WARN', sprintf('mixed algorithms: %s', strjoin(cellstr(algs), ', ')));
end

% Coverage of the subject x speed x leg grid, so a crashed fit shows up as a
% hole rather than as a silently smaller table.
for s = unique(T.subject)'
    sel = T.subject == s;
    fprintf('    subject %d: %d fits, speeds [%s], legs [%s], n = %d, %d samples x %d trials\n', ...
        s, sum(sel), num2str(unique(T.speed(sel))'), num2str(unique(T.leg(sel))'), ...
        T.n(find(sel, 1)), T.nsample(find(sel, 1)), T.ntrial(find(sel, 1)));
end

% ---------------------------------------------------------------------
fprintf('\n--- 2. CONVERGENCE --------------------------------------------------\n');

% fmincon exit flags: 1 = first-order optimality met, 2 = step below
% StepTolerance, 3 = objective change below FunctionTolerance, 0 = hit an
% evaluation cap, negatives = failure.  Only 1 is convergence in the usual
% sense; 2 means the line search stalled, which on this problem has meant a
% noisy gradient rather than a real minimum.
uef = unique(T.exitflag);
for k = 1:numel(uef)
    fprintf('    exitflag %+d  x%-3d  %s\n', uef(k), sum(T.exitflag == uef(k)), ...
        exitflag_label(uef(k)));
end

if all(T.exitflag == 1)
    verdict('PASS', 'every fit met the first-order optimality tolerance');
elseif any(T.exitflag < 0)
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf('%d fit(s) exited with a negative flag.', sum(T.exitflag < 0)));
elseif all(T.exitflag == 2)
    n_warn = n_warn + 1;
    verdict('WARN', ['every fit exited on flag 2 (step below StepTolerance).  fmincon ' ...
        'stopped moving, not because the gradient vanished.  Check firstorderopt.']);
else
    n_warn = n_warn + 1;
    verdict('WARN', 'not every fit reached flag 1.');
end

fprintf('    iterations    min %5d  median %6.0f  max %5d\n', ...
    min(T.iterations), median(T.iterations), max(T.iterations));
fprintf('    funcCount     min %5d  median %6.0f  max %5d\n', ...
    min(T.funcCount), median(T.funcCount), max(T.funcCount));
fprintf('    firstorderopt min %9.2e  median %9.2e  max %9.2e\n', ...
    min(T.firstorderopt), median(T.firstorderopt), max(T.firstorderopt));

if median(T.firstorderopt) > 1
    n_warn = n_warn + 1;
    verdict('WARN', sprintf(['median firstorderopt is %.2g -- the gradient is still ' ...
        'large at the reported solution, so these are stalls, not minima.'], ...
        median(T.firstorderopt)));
end

if all(isnan(T.time_s))
    verdict('INFO', ['wall-clock time was not saved (no t_fit_s).  Using funcCount ' ...
        'as the cost proxy; add tic/toc to main.m for future runs.']);
else
    fprintf('    time (s)      min %8.1f  median %8.1f  max %8.1f  total %8.1f (%.2f h)\n', ...
        min(T.time_s), median(T.time_s), max(T.time_s), ...
        sum(T.time_s), sum(T.time_s)/3600);
end

% ---------------------------------------------------------------------
fprintf('\n--- 3. SOLUTION QUALITY ---------------------------------------------\n');

fprintf('    objective E   min %8.3f  median %8.3f  max %8.3f\n', ...
    min(T.fval), median(T.fval), max(T.fval));
fprintf('    per-trial RMSE spread (max-min) within a fit: median %.3f, worst %.3f\n', ...
    median(T.rmse_max - T.rmse_min), max(T.rmse_max - T.rmse_min));

% The saved Q,l must reproduce the objective fmincon reported.  A mismatch
% means the reconstruction path (theta_to_Ql with the right cond_max) is
% wrong, which is exactly the bug that a changed eps_shift would cause.
tol_fval = 1e-6 * max(1, median(T.fval));
if all(T.fval_resid < tol_fval)
    verdict('PASS', sprintf('saved (Q,l) reproduces fval in every fit (max |diff| %.2e)', ...
        max(T.fval_resid)));
else
    n_fail = n_fail + 1;
    bad = find(T.fval_resid >= tol_fval);
    verdict('FAIL', sprintf(['re-solving the QP with the saved (Q,l) does not reproduce ' ...
        'fval in %d of %d fits (max |diff| %.2e).'], numel(bad), height(T), max(T.fval_resid)));
    for k = bad'
        fprintf('           subject %d speed %d leg %d:  fmincon %.3f  ->  re-solve %.3f  (cond(Q) = %.2e)\n', ...
            T.subject(k), T.speed(k), T.leg(k), T.fval(k), T.rmse_pooled(k), T.cond_Q(k));
    end
    fprintf(['           The same Q on the same data gives different forces, so the QP\n' ...
             '           minimiser is not numerically unique.  Suspect ill-conditioned Q,\n' ...
             '           loose OSQP tolerances, or both.\n']);
end

fprintf('    cond(Q)       min %9.3e  median %9.3e  max %9.3e\n', ...
    min(T.cond_Q), median(T.cond_Q), max(T.cond_Q));
fprintf('    lambda_min    min %9.3e  median %9.3e  max %9.3e\n', ...
    min(T.lambda_min), median(T.lambda_min), max(T.lambda_min));
fprintf('    trace(Q)/n    min %9.3e  median %9.3e  max %9.3e\n', ...
    min(T.trace_Q ./ T.n), median(T.trace_Q ./ T.n), max(T.trace_Q ./ T.n));
fprintf('    ||l||         min %9.3e  median %9.3e  max %9.3e\n', ...
    min(T.norm_l), median(T.norm_l), max(T.norm_l));

% Q must be positive definite for the DOC minimiser to be unique.  A
% non-positive lambda_min means it is not, and the fit rests on nothing.
n_bad_psd = sum(T.lambda_min <= 0);
if n_bad_psd == 0
    verdict('PASS', 'Q is positive definite in every fit');
else
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf(['lambda_min <= 0 in %d fit(s) -- Q is singular or ' ...
        'indefinite, so the QP minimiser is not unique and the fit is not identified.'], ...
        n_bad_psd));
end

% The bound the constraint is supposed to buy.  Checked against each file's
% own cond_max where recorded, and against the 1e4 default otherwise.
cmax = T.cond_max;
cmax(isnan(cmax)) = 1e4;
n_over = sum(T.cond_Q > cmax .* (1 + 1e-6));
if n_over == 0
    verdict('PASS', 'cond(Q) within cond_max in every fit');
else
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf('cond(Q) exceeds cond_max in %d of %d fits (worst %.3e).', ...
        n_over, height(T), max(T.cond_Q)));
end

if any(T.constrained)
    sel = T.constrained;
    n_off = sum(T.trace_resid(sel) > 1e-4 * T.n(sel));
    if n_off == 0
        verdict('PASS', sprintf('trace(Q) = n holds (max residual %.2e)', ...
            max(T.trace_resid(sel))));
    else
        n_warn = n_warn + 1;
        verdict('WARN', sprintf(['trace(Q) is off by more than 1e-4*n in %d fit(s) ' ...
            '(max %.2e) -- fmincon stopped before restoring feasibility.'], ...
            n_off, max(T.trace_resid(sel))));
    end
else
    % Without the equality the objective is invariant along (Q,l) -> (cQ,cl),
    % so trace(Q)/n measures how far the iterate slid along the gauge instead
    % of fitting.  It is a diagnostic, not a failure on its own.
    verdict('INFO', sprintf(['unconstrained run: trace(Q)/n ranges %.3g to %.3g, ' ...
        'i.e. the iterate drifted up to %.3gx along the (Q,l) -> (cQ,cl) gauge.'], ...
        min(T.trace_Q ./ T.n), max(T.trace_Q ./ T.n), max(T.trace_Q ./ T.n)));
end

% =====================================================================
% Per-fit table
% =====================================================================
fprintf('\n--- PER-FIT ---------------------------------------------------------\n');
fprintf('%-4s %-5s %-3s %4s %5s %6s %9s %8s %8s %10s %10s %9s\n', ...
    'subj', 'speed', 'leg', 'ef', 'iter', 'nfev', 'E', 'rmse_sd', 'time_s', ...
    'cond(Q)', 'lam_min', 'trace/n');
for k = 1:height(T)
    fprintf('%-4d %-5d %-3d %4d %5d %6d %9.4f %8.3f %8.1f %10.3e %10.3e %9.3e\n', ...
        T.subject(k), T.speed(k), T.leg(k), T.exitflag(k), T.iterations(k), ...
        T.funcCount(k), T.fval(k), T.rmse_std(k), T.time_s(k), ...
        T.cond_Q(k), T.lambda_min(k), T.trace_Q(k) / T.n(k));
end

fprintf('\n--- VERDICT ---------------------------------------------------------\n');
if n_fail > 0
    fprintf('  %d FAIL, %d WARN.  These fits should NOT be used.\n', n_fail, n_warn);
elseif n_warn > 0
    fprintf('  0 FAIL, %d WARN.  Usable, but read the warnings above.\n', n_warn);
else
    fprintf('  All checks passed.\n');
end
fprintf('=====================================================================\n\n');

% =====================================================================
% Optional outputs
% =====================================================================
if want_csv
    csv_path = fullfile(results_dir, 'summary.csv');
    writetable(T, csv_path);
    fprintf('  Wrote %s\n', csv_path);
end

if want_plot
    fig = summary_figure(T);
    png_path = fullfile(results_dir, 'summary.png');
    exportgraphics(fig, png_path, 'Resolution', 150);
    close(fig);
    fprintf('  Wrote %s\n', png_path);
end

end

% =========================================================================
% Helpers
% =========================================================================

function v = get_or(S, name, default)
%GET_OR  Field lookup that tolerates files written by older code.
if isstruct(S) && isfield(S, name)
    v = S.(name);
else
    v = default;
end
end

function verdict(tag, msg)
%VERDICT  One check, one line.  PASS / WARN / FAIL / INFO in the left column
%so the report can be skimmed for the word FAIL.
fprintf('  [%-4s] %s\n', tag, msg);
end

function s = exitflag_label(ef)
switch ef
    case  1, s = 'first-order optimality within tolerance';
    case  2, s = 'step smaller than StepTolerance';
    case  3, s = 'objective change smaller than FunctionTolerance';
    case  0, s = 'hit MaxIterations or MaxFunctionEvaluations';
    case -1, s = 'stopped by an output or plot function';
    case -2, s = 'no feasible point found / converged to an infeasible point';
    case -3, s = 'objective below ObjectiveLimit';
    otherwise, s = 'unrecognised exit flag';
end
end

function fig = summary_figure(T)
fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0 0 0.9 0.7]);
lbl = arrayfun(@(a, b, c) sprintf('s%d/%d/%d', a, b, c), ...
    T.subject, T.speed, T.leg, 'UniformOutput', false);
x = 1:height(T);

subplot(2, 2, 1);
bar(x, T.fval); grid on;
set(gca, 'XTick', x, 'XTickLabel', lbl, 'XTickLabelRotation', 90);
ylabel('objective E (RMSE)'); title('Fit quality');

subplot(2, 2, 2);
semilogy(x, max(T.cond_Q, 1), 'o-'); grid on; hold on;
cmax = T.cond_max; cmax(isnan(cmax)) = 1e4;
semilogy(x, cmax, 'r--', 'LineWidth', 1.2);
set(gca, 'XTick', x, 'XTickLabel', lbl, 'XTickLabelRotation', 90);
ylabel('cond(Q)'); title('Conditioning'); legend('cond(Q)', 'cond_{max}', 'Location', 'best');

subplot(2, 2, 3);
bar(x, T.iterations); grid on;
set(gca, 'XTick', x, 'XTickLabel', lbl, 'XTickLabelRotation', 90);
ylabel('fmincon iterations'); title('Optimisation effort');

subplot(2, 2, 4);
errorbar(x, T.rmse_mean, T.rmse_mean - T.rmse_min, T.rmse_max - T.rmse_mean, 'o'); grid on;
set(gca, 'XTick', x, 'XTickLabel', lbl, 'XTickLabelRotation', 90);
ylabel('per-trial RMSE'); title('Within-fit spread (min/mean/max)');
end
