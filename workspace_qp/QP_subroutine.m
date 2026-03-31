function [Fout, varargout] = QP_subroutine(Q, l, data, vars, model, sample_list, trial_list, speed_list, leg_list, varargin)
%QP_SUBROUTINE  Solve the parametric QP over all requested conditions.
%
%   Fout = QP_SUBROUTINE(Q, l, data, vars, model, ...)
%   [Fout, dFout] = QP_SUBROUTINE(Q, l, data, vars, model, ...)
%
%   When called with two outputs, also returns dFout: the sensitivity of
%   each optimal force vector f* w.r.t. the IO parameter vector
%   theta = [q; l], where q holds the lower-triangular entries of the
%   Cholesky factor L of Q (Q = L*L').
%
%   Outputs:
%     Fout  - [n x 1 x nsamples x ntrials x nspeeds x nlegs]
%     dFout - [n x ntheta x nsamples x ntrials x nspeeds x nlegs]  (optional)
%             ntheta = n*(n+1)/2 + n

sens_flag = nargout >= 2;

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
    if ~isempty(varargin)
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

                sol   = model.solve();
                f_opt = sol.value(vars.variables.f);

                if sens_flag
                    lam_g  = sol.value(model.lam_g);
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

end
