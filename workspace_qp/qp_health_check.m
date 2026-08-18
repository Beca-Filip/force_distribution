function report = qp_health_check(data_path, varargin)
%QP_HEALTH_CHECK  Diagnose an IO fit in a minute instead of a day.
%
%   report = QP_HEALTH_CHECK(data_path)
%   report = QP_HEALTH_CHECK(data_path, 'Theta', theta, ...)
%   report = QP_HEALTH_CHECK(data_path, 'Checkpoint', 'main/checkpoint.mat')
%
%   Everything that has gone wrong with these fits so far was decidable from
%   a dozen QPs at the right theta, and was instead discovered after hours of
%   remote-desktop batch: the CasADi abort on subject 4 speed 1 leg 2, and
%   the first-order optimality that oscillates while the objective falls.
%   This runs the same checks on a small subset and prints a verdict.
%
%   It is designed to be pointed at a RUNNING fit.  A fit writes
%   checkpoint.mat every iteration (see QP_LOG), so from a second MATLAB
%
%       qp_health_check('../Optimization Model Data/Patient4.mat', ...
%                       'Checkpoint', 'main/checkpoint.mat', ...
%                       'SpeedList', 1, 'LegList', 2);
%
%   answers "is this run worth continuing" without touching the run.
%
%
%   WHAT IT CHECKS
%   --------------
%   1. CONDITIONING.  cond(Sigma), cond(W) and cond(H) for H = W'*Q*W, the
%      matrix the QP solver actually sees.  The conditioning machinery bounds
%      cond(Q) <= cond_max; it does not bound cond(H), and the whitening
%      contributes cond(W)^2 on top -- 3.7e6 for patient 4, 1.8e7 for
%      patient 5, before Q does anything at all.  This block reports the
%      product and the ceiling it can reach.
%
%   2. SOLVER.  Every QP in the subset through QP_SOLVE, counting how many
%      needed a relaxed tolerance and how many did not solve at all.  A run
%      that is about to abort shows up here as ladder rungs being used.
%
%   3. GRADIENT.  Analytic dE/dtheta against central finite differences along
%      random directions, then per QP.  On subject 4 at theta0 the two agree
%      to 1e-6; after 120 fmincon iterations they agree only to 1e-2, and the
%      per-QP breakdown shows why -- a few QPs sit at an active-set kink,
%      where f*(theta) has no two-sided derivative at all.
%
%   4. DEGENERACY.  How many QPs are at such a kink, the smallest active
%      multiplier, rcond of the KKT matrix, and how often the minimum-norm
%      fallback fires.  This is the census that explains item 3.
%
%
%   NAME-VALUE ARGUMENTS
%     Theta        [], parameter vector to test at.  Default: theta0, i.e.
%                  Q = I with a seeded random l -- the point every fit starts
%                  from.
%     Checkpoint   '', path to a checkpoint.mat; its theta overrides Theta.
%     SampleList   1:20:101      subset, kept small on purpose
%     TrialList    1:2
%     SpeedList    1
%     LegList      1
%     CondMax      1e4
%     NDirections  3, random directions for the directional derivative check
%     Seed         0, for l0 and for the random directions
%     SolverOpt    [], solver option struct; default is main.m's tight OSQP
%
%   The returned struct carries every number printed, so this is usable as a
%   regression check as well as a console tool.
%
%   See also QP_SOLVE, QP_LOG, KKT_SENSITIVITY, SUMMARIZE_MAIN_RESULTS.

p = inputParser;
p.addParameter('Theta',       []);
p.addParameter('Checkpoint',  '');
p.addParameter('SampleList',  1:20:101);
p.addParameter('TrialList',   1:2);
p.addParameter('SpeedList',   1);
p.addParameter('LegList',     1);
p.addParameter('CondMax',     1e4);
p.addParameter('NDirections', 3);
p.addParameter('Seed',        0);
p.addParameter('SolverOpt',   []);
p.parse(varargin{:});
opt = p.Results;

% =========================================================================
% Setup
% =========================================================================
S    = load(data_path);
data = prepare_qp_data(S.data);
n    = size(data.f, 1);
[model, vars] = form_casadi_qp_model(n);

sol_opt = opt.SolverOpt;
if isempty(sol_opt)
    sol_opt = struct();
    % See main.m: this suppresses CasADi's dump of every solver input on a
    % failed QP and surfaces OSQP's return status instead.
    sol_opt.error_on_fail = false;
    sol_opt.osqp.verbose  = 0;
    sol_opt.osqp.eps_abs  = 1e-10;
    sol_opt.osqp.eps_rel  = 1e-10;
    sol_opt.osqp.max_iter = 200000;
end
model.solver('osqp', sol_opt);

cond_max  = opt.CondMax;
eps_shift = n / cond_max;
nq        = n * (n + 1) / 2;

theta  = opt.Theta;
origin = 'supplied theta';
if ~isempty(opt.Checkpoint)
    ck     = load(opt.Checkpoint);
    theta  = ck.theta;
    origin = sprintf('%s (iteration %d, fval %.6g)', opt.Checkpoint, ...
        get_or(ck, 'iteration', -1), get_or(ck, 'fval', NaN));
end
if isempty(theta)
    mask = tril(true(n));
    [r, c] = find(mask);
    q0 = zeros(nq, 1);
    q0(r == c) = sqrt(1 - eps_shift);
    rng(opt.Seed, 'twister');
    theta  = [q0; randn(n, 1)];
    origin = 'theta0 (Q = I, seeded random l)';
end
theta = reshape(theta, [], 1);

sample_list = opt.SampleList;
trial_list  = opt.TrialList;
speed_list  = opt.SpeedList;
leg_list    = opt.LegList;
n_qp = numel(sample_list) * numel(trial_list) * numel(speed_list) * numel(leg_list);

fprintf('\n=========================================================================\n');
fprintf(' QP HEALTH CHECK\n');
fprintf('   data      %s  (n = %d)\n', data_path, n);
fprintf('   theta     %s\n', origin);
fprintf('   subset    %d samples x %d trials x %d speeds x %d legs = %d QPs\n', ...
    numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list), n_qp);
fprintf('=========================================================================\n');

report = struct();
report.n_qp   = n_qp;
report.origin = origin;
n_fail = 0;
n_warn = 0;

% =========================================================================
% 1. Conditioning
% =========================================================================
fprintf('\n---- 1. CONDITIONING ----------------------------------------------------\n');

[Q, l, L] = theta_to_Ql(theta, n, cond_max);
W  = data.F_invcov;
H  = (W' * Q * W);
H  = (H + H') / 2;

eQ = sort(eig((Q + Q') / 2));
report.cond_Q     = eQ(end) / eQ(1);
report.trace_Q    = trace(Q);
report.cond_W     = cond(W);
report.cond_H     = cond(H);
report.cond_H_max = cond(W)^2 * cond_max;

fprintf('  cond(Q)        %.3e   (ceiling %.1e)\n', report.cond_Q, cond_max);
fprintf('  trace(Q)       %.4f     (should be exactly %d)\n', report.trace_Q, n);
fprintf('  cond(W)        %.3e   whitening alone\n', report.cond_W);
fprintf('  cond(H)        %.3e   H = W''*Q*W, what OSQP actually factorises\n', report.cond_H);
fprintf('  cond(H) worst  %.3e   = cond(W)^2 * cond_max, reachable within the constraint\n', ...
    report.cond_H_max);

% fmincon's interior-point algorithm satisfies a nonlinear equality only in
% the limit, so a checkpoint taken mid-run is legitimately off the trace
% constraint by whatever the constraint violation column says -- typically
% 1e-3 after a hundred iterations.  Only a gross violation says the theta was
% produced without the constraint at all; anything smaller is a note, and the
% cond(Q) test below is the one that carries the guarantee.
trace_violation = abs(report.trace_Q - n) / n;
report.trace_violation = trace_violation;

if trace_violation > 1e-2
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf(['trace(Q) = %.4f, not %d.  Too far off to be an ' ...
        'unconverged iterate: this theta was probably produced without the trace ' ...
        'constraint, in which case the cond(Q) bound does not hold for it.'], ...
        report.trace_Q, n));
elseif trace_violation > 1e-6
    n_warn = n_warn + 1;
    verdict('WARN', sprintf(['trace(Q) = %.4f against %d, a relative violation of ' ...
        '%.1e.  Expected for a checkpoint from a running fit; not expected for a ' ...
        'saved final theta.'], report.trace_Q, n, trace_violation));
elseif report.cond_Q > cond_max * (1 + 1e-6)
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf('cond(Q) = %.3e exceeds the ceiling %.1e.', ...
        report.cond_Q, cond_max));
else
    verdict('PASS', 'Q is trace-normalised and inside its condition ceiling.');
end

if report.cond_H > 1e8
    n_warn = n_warn + 1;
    verdict('WARN', sprintf(['cond(H) = %.2e.  OSQP is a first-order method; past ' ...
        'about 1e8 it stops reaching eps 1e-10 in any affordable iteration count, ' ...
        'and the KKT solve behind the gradient loses the same digits.'], report.cond_H));
else
    verdict('PASS', sprintf('cond(H) = %.2e, within what OSQP handles at eps 1e-10.', ...
        report.cond_H));
end

verdict('INFO', ['the conditioning guarantee bounds cond(Q), not cond(H).  ' ...
    sprintf('cond(W)^2 = %.2e is spent before Q contributes anything.', report.cond_W^2)]);

% =========================================================================
% 2. Solver census
% =========================================================================
fprintf('\n---- 2. QP SOLVER -------------------------------------------------------\n');

rungs   = zeros(n_qp, 1);
t_qp    = zeros(n_qp, 1);
solved  = true(n_qp, 1);
f_all   = zeros(n, n_qp);
lam_all = cell(n_qp, 1);
idx     = 0;
first_failure = '';

for trial = trial_list
for speed = speed_list
for leg   = leg_list
for k     = sample_list
    idx = idx + 1;
    model = set_qp_parameters(data, vars, model, k, trial, speed, leg);
    model = set_qp_weights(Q, l, vars, model);
    model = set_qp_normalization(data, vars, model);
    try
        [f_k, lam_k, sinfo] = qp_solve(model, vars, sol_opt);
        f_all(:, idx)  = f_k;
        lam_all{idx}   = lam_k;
        rungs(idx)     = sinfo.rung;
        t_qp(idx)      = sinfo.t_s;
    catch ME
        solved(idx) = false;
        if isempty(first_failure)
            first_failure = sprintf('sample %d trial %d speed %d leg %d: %s', ...
                k, trial, speed, leg, ME.message);
        end
    end
end
end
end
end

report.n_unsolved = sum(~solved);
report.n_degraded = sum(rungs > 0);
report.t_qp_median = median(t_qp(solved));
report.t_qp_max    = max(t_qp(solved));

fprintf('  solved            %d / %d\n', sum(solved), n_qp);
fprintf('  needed a retry    %d  (rungs used: %s)\n', report.n_degraded, ...
    mat2str(unique(rungs(rungs > 0))'));
fprintf('  time per QP       median %.2f ms, worst %.2f ms\n', ...
    1e3 * report.t_qp_median, 1e3 * report.t_qp_max);

if report.n_unsolved > 0
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf(['%d of %d QPs did not solve at any tolerance.  ' ...
        'A batch run would have aborted here.\n           first: %s'], ...
        report.n_unsolved, n_qp, first_failure));
elseif report.n_degraded > 0
    n_warn = n_warn + 1;
    verdict('WARN', sprintf(['%d of %d QPs needed a relaxed tolerance.  They are ' ...
        'solved, but their KKT sensitivity is correspondingly less accurate, and ' ...
        'the count rises as cond(H) grows -- this is the run heading for the abort.'], ...
        report.n_degraded, n_qp));
else
    verdict('PASS', 'every QP solved at the requested tolerance, first attempt.');
end

if ~all(solved)
    fprintf('\n(stopping: the gradient checks below need every QP to solve)\n');
    report.n_fail = n_fail;
    report.n_warn = n_warn;
    return;
end

% =========================================================================
% 3. Gradient
% =========================================================================
fprintf('\n---- 3. GRADIENT --------------------------------------------------------\n');

fun = @(th) QP_IO_inner_loop(th, data, vars, model, ...
    sample_list, trial_list, speed_list, leg_list, cond_max, sol_opt);

[E0, g0] = fun(theta);
fprintf('  E = %.8g,  ||grad|| = %.4e\n', E0, norm(g0));

rng(opt.Seed + 1, 'twister');
h_list   = [1e-4, 1e-5, 1e-6];
rel_err  = zeros(opt.NDirections, numel(h_list));
for di = 1:opt.NDirections
    d  = randn(numel(theta), 1);
    d  = d / norm(d);
    gd = g0' * d;
    fprintf('  direction %d:  analytic %+.6e\n', di, gd);
    for hi = 1:numel(h_list)
        h  = h_list(hi);
        fd = (fun(theta + h * d) - fun(theta - h * d)) / (2 * h);
        rel_err(di, hi) = abs(fd - gd) / max(abs(gd), 1e-12);
        fprintf('                h = %.0e   FD %+.6e   rel.err %.2e\n', h, fd, rel_err(di, hi));
    end
end

% The best rel.err over h is the honest figure: a large h is truncation and a
% small h is round-off, and neither is a statement about the formula.
report.grad_rel_err = min(rel_err(:, :), [], 2);
worst = max(report.grad_rel_err);
fprintf('  worst directional relative error (best h per direction): %.2e\n', worst);

if worst < 1e-4
    verdict('PASS', sprintf(['analytic and finite-difference directional derivatives ' ...
        'agree to %.1e.  The gradient is right here.'], worst));
elseif worst < 1e-2
    n_warn = n_warn + 1;
    verdict('WARN', sprintf(['gradient and finite differences disagree by %.1e.  ' ...
        'That is too large to be round-off and too small to be a wrong formula; ' ...
        'see the per-QP breakdown below.'], worst));
else
    n_fail = n_fail + 1;
    verdict('FAIL', sprintf('gradient and finite differences disagree by %.1e.', worst));
end

% =========================================================================
% 4. Degeneracy census -- why item 3 came out the way it did
% =========================================================================
fprintf('\n---- 4. ACTIVE-SET DEGENERACY ------------------------------------------\n');

n_act     = zeros(n_qp, 1);
rcond_K   = zeros(n_qp, 1);
min_mult  = zeros(n_qp, 1);
min_gap   = zeros(n_qp, 1);
degen     = false(n_qp, 1);
fallback  = false(n_qp, 1);
qp_rel    = nan(n_qp, 1);
fd_norm   = nan(n_qp, 1);

% One shared direction, so the per-QP errors are comparable with each other.
rng(opt.Seed + 2, 'twister');
d = randn(numel(theta), 1);
d = d / norm(d);
h = 1e-6;
[Qp, lp] = theta_to_Ql(theta + h * d, n, cond_max);
[Qm, lm] = theta_to_Ql(theta - h * d, n, cond_max);

idx = 0;
for trial = trial_list
for speed = speed_list
for leg   = leg_list
for k     = sample_list
    idx = idx + 1;

    [df, kinfo] = kkt_sensitivity(f_all(:, idx), lam_all{idx}, Q, L, data, k, trial, speed, leg);

    n_act(idx)    = kinfo.n_act;
    rcond_K(idx)  = kinfo.rcond_K;
    min_mult(idx) = kinfo.min_mult;
    min_gap(idx)  = kinfo.min_gap;
    degen(idx)    = kinfo.degenerate;
    fallback(idx) = kinfo.fallback;

    fp = solve_at(Qp, lp, data, vars, model, sol_opt, k, trial, speed, leg);
    fm = solve_at(Qm, lm, data, vars, model, sol_opt, k, trial, speed, leg);
    fd = (fp - fm) / (2 * h);

    % Scaled by the derivative's own size OR by a small fraction of the
    % force scale, whichever is larger.  A QP pinned at a vertex has
    % df/dtheta ~ 0 along most directions; dividing by that alone turns a
    % numerically irrelevant discrepancy into a headline 50% error, which is
    % how the first version of this check managed to condemn theta0.
    fd_norm(idx) = norm(fd);
    qp_rel(idx)  = norm(df * d - fd) / max([norm(fd), 1e-4 * norm(f_all(:, idx)), 1e-12]);
end
end
end
end

report.n_act      = n_act;
report.rcond_K    = rcond_K;
report.min_mult   = min_mult;
report.degenerate = degen;
report.qp_rel_err = qp_rel;
report.fd_norm    = fd_norm;

fprintf('  active bounds per QP    median %g, max %g\n', median(n_act), max(n_act));
fprintf('  rcond(KKT matrix)       median %.2e, worst %.2e\n', median(rcond_K), min(rcond_K));
fprintf('  smallest active mult.   median %.2e, worst %.2e\n', median(min_mult), min(min_mult));
fprintf('  minimum-norm fallback   %d / %d QPs\n', sum(fallback), n_qp);
fprintf('  per-QP sensitivity error vs FD:  median %.2e, worst %.2e\n', ...
    median(qp_rel), max(qp_rel));

bad = find(qp_rel > 1e-3);
if ~isempty(bad)
    fprintf('  QPs above 1e-3:\n');
    for bi = bad(:)'
        fprintf('     #%-3d  rel.err %.2e   n_act %2d   min mult %.2e   rcond(K) %.2e   %s\n', ...
            bi, qp_rel(bi), n_act(bi), min_mult(bi), rcond_K(bi), ...
            ternary_str(degen(bi), 'DEGENERATE', ''));
    end
end

if isempty(bad)
    verdict('PASS', 'every QP''s sensitivity matches finite differences.');
else
    n_warn = n_warn + 1;
    verdict('WARN', sprintf(['%d of %d QPs have a sensitivity error above 1e-3.  ' ...
        'At those points an active bound is held with a near-zero multiplier, so ' ...
        'f*(theta) has a kink and no two-sided derivative exists.  The formula is ' ...
        'not wrong; the function is not differentiable there.'], numel(bad), n_qp));
    verdict('INFO', ['pooled over the full 1010 QPs of a real fit, some QP is at a ' ...
        'kink at essentially every theta.  fmincon assumes a smooth objective, so ' ...
        'the expected signature is exactly what the B5 run showed: the objective ' ...
        'falls monotonically, the step shrinks, and first-order optimality ' ...
        'oscillates instead of converging.  Treat a run that stops on step ' ...
        'tolerance as finished, not as broken.']);
end

% =========================================================================
fprintf('\n=========================================================================\n');
fprintf(' %d failures, %d warnings\n', n_fail, n_warn);
fprintf('=========================================================================\n');
report.n_fail = n_fail;
report.n_warn = n_warn;

end


% =========================================================================
function f = solve_at(Q, l, data, vars, model, sol_opt, k, trial, speed, leg)
model = set_qp_parameters(data, vars, model, k, trial, speed, leg);
model = set_qp_weights(Q, l, vars, model);
model = set_qp_normalization(data, vars, model);
f = qp_solve(model, vars, sol_opt);
end


function verdict(tag, msg)
fprintf('  [%-4s] %s\n', tag, msg);
end


function s = ternary_str(cond_, a, b)
if cond_
    s = a;
else
    s = b;
end
end


function v = get_or(S, name, default)
if isstruct(S) && isfield(S, name)
    v = S.(name);
else
    v = default;
end
end
