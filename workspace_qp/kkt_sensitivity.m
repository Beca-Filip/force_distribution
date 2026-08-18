function [df_dtheta, info] = kkt_sensitivity(f_opt, lam_g, Q, L, data, k, trial, speed, leg)
%KKT_SENSITIVITY  Sensitivity of the QP solution w.r.t. the cost parameters.
%
%   df_dtheta        = KKT_SENSITIVITY(f_opt, lam_g, Q, L, data, k, trial, speed, leg)
%   [df_dtheta, info] = KKT_SENSITIVITY(...)
%
%   Uses the implicit function theorem applied to the KKT conditions of the
%   parametric QP.  At the solution (f*, nu*, lambda*), differentiating the
%   KKT stationarity condition w.r.t. theta = [q; l] gives the linear system
%
%       K * [df*/dtheta; dnu/dtheta; dlambda_act/dtheta] = -dR/dtheta
%
%   where K is the (reduced) KKT matrix and R is the stationarity residual.
%   The first n rows of the solution give df*/dtheta.
%
%   The cost is  0.5*f_norm'*Q*f_norm + l'*f_norm,  f_norm = W*(f - mu),
%   where W = data.F_invcov and mu = data.f_mean.
%
%   Inputs:
%     f_opt   - [n x 1]          optimal primal solution
%     lam_g   - [(ne+2n) x 1]    dual variables from sol.value(model.lam_g):
%                                  lam_g(1:ne)       equality  (A*f - b = 0)
%                                  lam_g(ne+1:ne+n)  lower bound  (fmin-f <= 0)
%                                  lam_g(ne+n+1:end) upper bound  (f-fmax <= 0)
%     Q       - [n x n]          current cost weight matrix
%     L       - [n x n]          lower-triangular Cholesky factor  (Q = L*L')
%     data    - data struct with fields f_mean, F_invcov, A, b, fmin, fmax
%     k, trial, speed, leg       indices into the data arrays
%
%   Outputs:
%     df_dtheta - [n x ntheta]   Jacobian of f* w.r.t. theta = [q; l],
%                                ntheta = n*(n+1)/2 + n
%     info      - struct of diagnostics, for QP_HEALTH_CHECK and the run log:
%                   n_act        number of active bound constraints
%                   rcond_K      reciprocal condition estimate of the KKT matrix
%                   cond_H       condition number of the cost Hessian W'*Q*W
%                   min_mult     smallest strictly positive active multiplier
%                   min_gap      smallest strictly positive primal bound gap
%                   degenerate   true when min_mult or min_gap sits at the
%                                identification tolerance, i.e. the active set
%                                is about to change and df_dtheta is a one-sided
%                                derivative only
%                   fallback     true when the minimum-norm solve was used
%
%
%   ACCURACY, AND WHEN IT FAILS
%   ---------------------------
%   Measured against central finite differences on subject 4 (speed 1, leg 2),
%   this Jacobian is right to about 1e-6 relative on a QP with a clean active
%   set -- and only to 1e-2, occasionally 1e-1, on a QP whose active set is
%   degenerate (an active bound held with a multiplier near zero, or a bound
%   that enters or leaves under an infinitesimal step).  That is not a defect
%   of the formula: at such a point f*(theta) genuinely has a kink and no
%   two-sided derivative exists.  info.degenerate flags it.
%
%   The pooled objective sums over 1010 QPs per condition, so a handful of
%   them sit at a kink at essentially every theta.  The consequence is
%   visible in any long fmincon run: the objective decreases monotonically and
%   the step size shrinks, while first-order optimality oscillates instead of
%   converging.  See QP_HEALTH_CHECK, which measures all of this in a minute
%   rather than over a 22-hour batch.

n  = length(f_opt);
nq = n*(n+1)/2;
ntheta = nq + n;
ne = size(data.A, 1);

W  = data.F_invcov;
mu = data.f_mean;
f_norm_star = W * (f_opt - mu);   % [n x 1]

A_val    = squeeze(data.A   (:, :, k, trial, speed, leg));   % [ne x n]
fmin_val = squeeze(data.fmin(:, :, k, trial, speed, leg));   % [n x 1]
fmax_val = squeeze(data.fmax(:, :, k, trial, speed, leg));   % [n x 1]

% ---- Dual variables -----------------------------------------------------
nu           = lam_g(1:ne);
lambda_lower = lam_g(ne+1    : ne+n);
lambda_upper = lam_g(ne+n+1  : ne+2*n);

% ---- Active-set identification ------------------------------------------
% A constraint is active if the primal gap is small OR the dual is non-zero.
tol = 1e-6;
active_lower = (f_opt - fmin_val < tol) | (lambda_lower > tol);
active_upper = (fmax_val - f_opt < tol) | (lambda_upper > tol);
n_act = sum(active_lower) + sum(active_upper);

% ---- KKT matrix ---------------------------------------------------------
H  = W' * Q * W;                          % [n x n] cost Hessian w.r.t. f
In = eye(n);
C_act = [-In(active_lower, :);
          In(active_upper, :)];            % [n_act x n]

K = [H,       A_val',              C_act';
     A_val,   zeros(ne,   ne),     zeros(ne,   n_act);
     C_act,   zeros(n_act, ne),    zeros(n_act, n_act)];   % [(n+ne+n_act) x same]

% ---- RHS: -d(stationarity)/dtheta --------------------------------------
%
% stationarity = W'*(Q*f_norm* + l) + A'*nu - lambda_lower + lambda_upper
%              = W'*Q*W*(f*-mu) + W'*l + (fixed dual terms)
%
% Only Q (via q) and l vary with theta; dual terms are treated as fixed.
%
% For l_k:   d(stationarity)/dl_k = W'[:,k]
%   => RHS_l = W'                                        [n x n]
%
% For q_k = L_{ij}  (i >= j):
%   dQ/dL_{ij} = e_i*(L[:,j])' + L[:,j]*e_i'
%   d(W'*Q*W*(f*-mu))/dL_{ij}
%     = W' * dQ/dL_{ij} * W*(f*-mu)
%     = W' * [e_i*(L[:,j])'*f_norm* + L[:,j]*(e_i'*f_norm*)]
%     = (L[:,j]'*f_norm*) * W'[:,i]  +  f_norm*[i] * W'*L[:,j]
%
% Precompute B = W'*L  =>  B[:,j] = W'*L[:,j]

B = W' * L;   % [n x n]

mask = tril(true(n));
[rows_lt, cols_lt] = find(mask);   % column-major lower-tri index pairs

% Vectorised form of
%     RHS_q(:,ki) = (L(:,jj)'*f_norm_star) * W(ii,:)' + f_norm_star(ii) * B(:,jj)
% column by column.  alpha depends only on the COLUMN index jj and beta only
% on the ROW index ii, so both are lookups into precomputed n-vectors and the
% whole [n x nq] block is two indexed products.  Bit-identical to the loop,
% about twice as fast, and this runs once per QP per gradient evaluation.
alpha = L' * f_norm_star;                                   % [n x 1], alpha(j)
RHS_q = W(rows_lt, :)' .* alpha(cols_lt)' ...
      + B(:, cols_lt)  .* f_norm_star(rows_lt)';            % [n x nq]

% Assemble full RHS  [(n+ne+n_act) x ntheta]
% The equality and active-inequality rows have zero sensitivity
% (A*f=b and C_act*f=g_act do not depend on theta).
RHS = -[RHS_q,                    W';
        zeros(ne + n_act, ntheta)];

% ---- Solve and extract primal sensitivity ------------------------------
% K can be singular when the active set is degenerate (e.g. both lower and
% upper bounds active simultaneously).  lsqminnorm gives the minimum-norm
% least-squares solution and does not produce Inf/NaN in that case.
%
% The deficiency test is rcond, not rank.  rank() takes an SVD of K on every
% one of the 1010 QPs per gradient evaluation (0.14 ms against rcond's 0.017);
% the threshold below, 1e-14, is the same relative singular value that rank()
% uses by default at this size (max(size(K))*eps ~ 1.7e-14), so the branch is
% taken on the same matrices.
rcond_K = rcond(K);
use_fallback = ~isfinite(rcond_K) || rcond_K < 1e-14;
if ~use_fallback
    sensitivity = K \ RHS;
    % A solve that squeaks past the rcond test can still return non-finite
    % entries.  Catching it here keeps NaNs out of the gradient, where they
    % would surface much later as an uninformative fmincon failure.
    use_fallback = ~all(isfinite(sensitivity(:)));
end
if use_fallback
    sensitivity = lsqminnorm(K, RHS);
end
df_dtheta   = sensitivity(1:n, :);   % [n x ntheta]

% ---- Diagnostics --------------------------------------------------------
if nargout > 1
    mult = [lambda_lower(active_lower); lambda_upper(active_upper)];
    gaps = [f_opt - fmin_val; fmax_val - f_opt];

    info = struct();
    info.n_act      = n_act;
    info.rcond_K    = rcond_K;
    info.cond_H     = cond(H);
    info.min_mult   = min_positive(mult, tol);
    info.min_gap    = min_positive(gaps, tol);
    info.fallback   = use_fallback;
    % "About to change" is judged on the same scale as the identification
    % tolerance: a multiplier or a gap within a factor of 1e3 of tol is close
    % enough that a fmincon step will move it across.
    info.degenerate = (info.min_mult < 1e3 * tol) || (info.min_gap < 1e3 * tol);
end

end


function m = min_positive(v, tol)
%MIN_POSITIVE  Smallest entry of v above tol, or Inf when there is none.
v = v(v > tol);
if isempty(v)
    m = Inf;
else
    m = min(v);
end
end
