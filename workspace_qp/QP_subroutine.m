function [Fout, varargout] = QP_subroutine(Q, l, data, vars, model, sample_list, trial_list, speed_list, leg_list, varargin)
%QP_SUBROUTINE  Solve the parametric QP over all requested conditions.
%
%   Fout            = QP_SUBROUTINE(Q, l, data, vars, model, ...)
%   [Fout, dFout]   = QP_SUBROUTINE(Q, l, data, vars, model, ...)
%   [..., solve_info] = QP_SUBROUTINE(...)
%
%   When called with two outputs, also returns dFout: the sensitivity of
%   each optimal force vector f* w.r.t. the IO parameter vector
%   theta = [q; l], where q holds the lower-triangular entries of the
%   Cholesky factor L of Q (Q = L*L').
%
%   Optional trailing arguments:
%     varargin{1}  L, the unshifted Cholesky factor, to avoid re-factorising
%                  Q at a trial point where it is PSD but not strictly PD
%     varargin{2}  the solver option struct the model was configured with.
%                  When given, every QP goes through QP_SOLVE, which retries
%                  on a tolerance ladder rather than letting a single hard QP
%                  raise a CasADi error and end the batch.  Without it the
%                  behaviour is the historical bare model.solve().
%
%   Outputs:
%     Fout       - [n x 1 x nsamples x ntrials x nspeeds x nlegs]
%     dFout      - [n x ntheta x nsamples x ntrials x nspeeds x nlegs]
%     solve_info - struct: n_degraded, worst_rung, t_solve_s.  Zeroed when no
%                  solver options were passed, since nothing was measured.

sens_flag = nargout >= 2;

if numel(varargin) >= 2 && ~isempty(varargin{2})
    sol_opt   = varargin{2};
    safe_mode = true;
else
    sol_opt   = [];
    safe_mode = false;
end

n_degraded = 0;
worst_rung = 0;
t_solve    = 0;

n        = size(vars.variables.f, 1);
nsamples = length(sample_list);
ntrials  = length(trial_list);
nspeeds  = length(speed_list);
nlegs    = length(leg_list);

Fout = zeros([n, 1, nsamples, ntrials, nspeeds, nlegs]);

if sens_flag
    ntheta = n*(n+1)/2 + n;
    dFout  = zeros([n, ntheta, nsamples, ntrials, nspeeds, nlegs]);
    % Use caller-supplied L when available (avoids re-factorising Q, which
    % can fail when Q is PSD but not strictly PD at a fmincon trial point).
    if numel(varargin) >= 1 && ~isempty(varargin{1})
        L = varargin{1};
    else
        L = chol(Q + 1e-12 * eye(n), 'lower');
    end
end

cntsamples = 1;
cnttrials  = 1;
cntspeeds  = 1;
cntlegs    = 1;

for trial = trial_list

    cntspeeds = 1;
    for speed = speed_list

        cntlegs = 1;
        for leg = leg_list

            cntsamples = 1;
            for k = sample_list

                model = set_qp_parameters(data, vars, model, k, trial, speed, leg);
                model = set_qp_weights(Q, l, vars, model);
                model = set_qp_normalization(data, vars, model);

                if safe_mode
                    % The context has to be attached here.  QP_SOLVE knows the
                    % solve failed but not which of the ten thousand QPs in a
                    % pooled evaluation it was, and that index is the whole
                    % difference between a reproducible bug and a lost night.
                    try
                        [f_opt, lam_g, sinfo] = qp_solve(model, vars, sol_opt);
                    catch ME
                        throw(MException('QP_subroutine:qpFailed', ...
                            ['QP failed at sample %d, trial %d, speed %d, leg %d ' ...
                             '(cond(Q) = %.3e, cond(W''*Q*W) = %.3e).\n%s'], ...
                            k, trial, speed, leg, cond(Q), ...
                            cond(data.F_invcov' * Q * data.F_invcov), ME.message));
                    end
                    n_degraded = n_degraded + double(sinfo.degraded);
                    worst_rung = max(worst_rung, sinfo.rung);
                    t_solve    = t_solve + sinfo.t_s;
                else
                    sol   = model.solve();
                    f_opt = sol.value(vars.variables.f);
                    if sens_flag
                        lam_g = sol.value(model.lam_g);
                    end
                end

                if sens_flag
                    df_dth = kkt_sensitivity(f_opt, lam_g, Q, L, data, k, trial, speed, leg);
                    dFout(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = df_dth;
                end

                Fout(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = f_opt;
                cntsamples = cntsamples + 1;
            end

            cntlegs = cntlegs + 1;
        end

        cntspeeds = cntspeeds + 1;
    end

    cnttrials = cnttrials + 1;
end

if sens_flag
    varargout{1} = dFout;
end

if nargout >= 3
    varargout{2} = struct('n_degraded', n_degraded, ...
                          'worst_rung', worst_rung, ...
                          't_solve_s',  t_solve);
end

end
