function [f_opt, lam_g, info] = qp_solve(model, vars, base_opt, ladder)
%QP_SOLVE  Solve the parametric QP once, without letting CasADi abort the run.
%
%   [f_opt, lam_g]       = QP_SOLVE(model, vars, base_opt)
%   [f_opt, lam_g, info] = QP_SOLVE(model, vars, base_opt, ladder)
%
%   The model must already carry its parameter values (SET_QP_PARAMETERS,
%   SET_QP_WEIGHTS, SET_QP_NORMALIZATION) and its solver must already be set
%   to base_opt.  This function only owns the solve.
%
%
%   WHY THIS EXISTS
%   ---------------
%   casadi.Opti.solve() THROWS when the underlying solver does not reach its
%   tolerance.  Inside a 22-hour batch that is fatal: one hard QP out of the
%   millions solved in a fit takes the whole run with it, and the error is a
%   SWIG:RuntimeError carrying a screenful of raw solver inputs rather than
%   anything that says which condition, sample or trial it was.
%
%   That is how the second B5 attempt died on subject 4, speed 1, leg 2.  The
%   mechanism is not a constraint problem -- the constraints are bounds plus
%   five linear equalities and satisfy LICQ throughout.  It is OSQP running
%   out of iterations, because the Hessian the QP actually sees is
%   H = W'*Q*W, and cond(H) grows as fmincon walks Q:
%
%       theta0     (Q = I)          cond(H) = 3.7e6   ~3 ms, solved easily
%       120 iters  (cond(Q) = 316)  cond(H) = 1.4e7   ~6 ms, needs more than
%                                                     4e3 OSQP iterations
%
%   and at the cond(Q) = 1e4 ceiling the trace constraint permits, cond(H)
%   reaches ~4e10, where OSQP does not reach eps_abs = 1e-10 within any
%   affordable iteration budget.  Note that cond(W)^2 = 3.7e6 is already
%   spent before Q contributes anything: the conditioning guarantee bounds
%   cond(Q), not cond(H), and the whitening sits in between.
%
%   So the failure is a tolerance failure, and the response is a tolerance
%   ladder.  Each rung re-solves the SAME model with a weaker request; a
%   retry after swapping the solver options was checked to reproduce a fresh
%   model's answer bit for bit.  A solve that needed a rung is reported in
%   info, never hidden -- a QP answered at eps 1e-6 instead of 1e-10 moves
%   the forces by about 3 parts in 1e3, which is tolerable for the objective
%   and is not tolerable for the KKT sensitivity that supplies the gradient.
%
%   opti.debug.value is deliberately not used as a last resort: after a conic
%   failure Opti has no iterate to hand back and raises on the attempt, so
%   there is nothing to salvage and nothing to be tempted by.
%
%   Inputs:
%     model     - the casadi.Opti object, parameters already set
%     vars      - the handle struct from FORM_CASADI_QP_MODEL
%     base_opt  - the solver option struct the model was configured with;
%                 restored before returning, whatever happened
%     ladder    - (optional) cell array of option structs to try in order
%                 after the base attempt fails.  Default: DEFAULT_LADDER
%                 below.  Pass {} to disable retrying and fail immediately.
%
%   Outputs:
%     f_opt  - [n x 1]        optimal forces
%     lam_g  - [(ne+2n) x 1]  constraint multipliers, as KKT_SENSITIVITY wants
%     info   - struct with
%                rung      0 when the base options succeeded, k for the kth
%                          ladder entry
%                degraded  true when rung > 0, i.e. this solution is less
%                          accurate than the run asked for
%                eps_used  the eps_abs actually used
%                status    solver return status string
%                t_s       wall time for the solve, retries included
%
%   Raises qp_solve:allRungsFailed, with every rung's status in the message,
%   when nothing works.  Callers log that and move on to the next condition
%   rather than losing the batch.
%
%   See also QP_SUBROUTINE, QP_HEALTH_CHECK, QP_LOG.

if nargin < 4
    ladder = default_ladder(base_opt);
end

t_start = tic;

info = struct('rung', 0, 'degraded', false, ...
              'eps_used', get_eps(base_opt), 'status', 'solved', 't_s', 0);

% ---- Base attempt -------------------------------------------------------
try
    sol   = model.solve();
    f_opt = sol.value(vars.variables.f);
    lam_g = sol.value(model.lam_g);
    info.t_s = toc(t_start);
    return;
catch ME
    base_reason = failure_reason(ME);
end

% ---- Ladder -------------------------------------------------------------
% base_opt is restored on every exit path, including the error one, so the
% caller's model is never left holding a relaxed tolerance that would then
% silently degrade every subsequent QP.
reasons = {sprintf('base (eps %.0e): %s', get_eps(base_opt), base_reason)};

for rung = 1:numel(ladder)
    try
        model.solver('osqp', ladder{rung});
        sol   = model.solve();
        f_opt = sol.value(vars.variables.f);
        lam_g = sol.value(model.lam_g);

        model.solver('osqp', base_opt);

        info.rung     = rung;
        info.degraded = true;
        info.eps_used = get_eps(ladder{rung});
        info.status   = 'solved';
        info.t_s      = toc(t_start);
        return;
    catch ME
        reasons{end+1} = sprintf('rung %d (eps %.0e, max_iter %g): %s', ...
            rung, get_eps(ladder{rung}), get_max_iter(ladder{rung}), ...
            failure_reason(ME));  %#ok<AGROW>
    end
end

model.solver('osqp', base_opt);

error('qp_solve:allRungsFailed', ...
      'the QP did not solve at any tolerance:\n           %s', ...
      strjoin(reasons, sprintf('\n           ')));

end


% =========================================================================
function ladder = default_ladder(base_opt)
%DEFAULT_LADDER  More iterations first, then less accuracy.
%
%   Ordered by what each rung costs the gradient.  Raising max_iter costs
%   time and nothing else, so it goes first; relaxing eps costs solution
%   accuracy, so it goes last and in two steps rather than one, to keep the
%   loosest rung rare and conspicuous in the log.
eps0 = get_eps(base_opt);
it0  = get_max_iter(base_opt);

ladder = { relax(base_opt, eps0,         10 * it0), ...
           relax(base_opt, 1e2 * eps0,   10 * it0), ...
           relax(base_opt, 1e4 * eps0,   10 * it0) };
end


function opt = relax(base_opt, eps_new, max_iter_new)
opt = base_opt;
if ~isfield(opt, 'osqp')
    opt.osqp = struct();
end
opt.osqp.eps_abs  = eps_new;
opt.osqp.eps_rel  = eps_new;
opt.osqp.max_iter = max_iter_new;
end


function e = get_eps(opt)
if isfield(opt, 'osqp') && isfield(opt.osqp, 'eps_abs')
    e = opt.osqp.eps_abs;
else
    e = 1e-3;                       % OSQP's own default
end
end


function m = get_max_iter(opt)
if isfield(opt, 'osqp') && isfield(opt.osqp, 'max_iter')
    m = opt.osqp.max_iter;
else
    m = 4000;                       % OSQP's own default
end
end


function r = failure_reason(ME)
%FAILURE_REASON  One line out of CasADi's multi-screen SWIG error.
%
%   Two shapes occur.  With error_on_fail left true (the default) the conic
%   interface raises first and the message says only "conic process failed";
%   with it false, Opti raises instead and the message carries the actual
%   OSQP status.  The status is the useful half, so prefer it when present.
tok = regexp(ME.message, 'return_status is ''([^'']*)''', 'tokens', 'once');
if ~isempty(tok)
    r = tok{1};
elseif contains(ME.message, 'conic process failed')
    r = 'conic process failed (status suppressed by error_on_fail)';
else
    lines = strsplit(ME.message, newline);
    r = strtrim(lines{end});
end
end
