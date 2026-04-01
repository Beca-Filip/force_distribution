function df_dtheta = kkt_sensitivity(f_opt, lam_g, Q, L, data, k, trial, speed, leg)
%KKT_SENSITIVITY  Sensitivity of the QP solution w.r.t. the cost parameters.
%
%   df_dtheta = KKT_SENSITIVITY(f_opt, lam_g, Q, L, data, k, trial, speed, leg)
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
%   Output:
%     df_dtheta - [n x ntheta]   Jacobian of f* w.r.t. theta = [q; l],
%                                ntheta = n*(n+1)/2 + n

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

RHS_q = zeros(n, nq);
for ki = 1:nq
    ii = rows_lt(ki);
    jj = cols_lt(ki);
    alpha = L(:, jj)' * f_norm_star;    % scalar: L[:,j]'*f_norm*
    beta  = f_norm_star(ii);             % scalar: f_norm*[i]
    RHS_q(:, ki) = alpha * W(ii, :)' + beta * B(:, jj);
end

% Assemble full RHS  [(n+ne+n_act) x ntheta]
% The equality and active-inequality rows have zero sensitivity
% (A*f=b and C_act*f=g_act do not depend on theta).
RHS = -[RHS_q,                    W';
        zeros(ne + n_act, ntheta)];

% ---- Solve and extract primal sensitivity ------------------------------
% K can be singular when the active set is degenerate (e.g. both lower and
% upper bounds active simultaneously).  lsqminnorm gives the minimum-norm
% least-squares solution and does not produce Inf/NaN in that case.
if rank(K) < size(K, 1)
    sensitivity = lsqminnorm(K, RHS);
else
    sensitivity = K \ RHS;
end
df_dtheta   = sensitivity(1:n, :);   % [n x ntheta]

end
