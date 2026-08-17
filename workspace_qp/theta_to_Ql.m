function [Q, l, L, eps_shift] = theta_to_Ql(theta, n, cond_max)
%THETA_TO_QL  Map the IO parameter vector to the QP cost weights (Q, l).
%
%   [Q, l, L, eps_shift] = THETA_TO_QL(theta, n)
%   [Q, l, L, eps_shift] = THETA_TO_QL(theta, n, cond_max)
%
%   theta = [q; l] with q the n*(n+1)/2 lower-triangular entries of L, and
%
%       Q = L*L' + eps_shift*I,     eps_shift = n / cond_max.
%
%   cond_max defaults to 1e4.  Pass Inf for the legacy unshifted Q = L*L'.
%
%   This is the ONLY place theta is turned into (Q, l).  Everything that needs
%   the mapping must come through here -- the objective, the search, and the
%   reconstruction of the saved result -- because a caller that rebuilt Q
%   without the shift would report a different matrix from the one that was
%   actually optimised.
%
%
%   WHY THE SHIFT AND THE TRACE CONSTRAINT
%   --------------------------------------
%   Two separate defects of the plain Q = L*L' parameterisation:
%
%   1. Scale gauge.  The QP minimiser is invariant under (Q,l) -> (cQ, cl) for
%      any c > 0: the cost scales by c and the constraints do not involve Q at
%      all.  So (Q,l) is identified only up to positive scale and the IO
%      objective has an exact flat direction.  fmincon spends its budget
%      sliding along it.
%
%   2. Conditioning.  Nothing stops L from becoming near-singular.  A Q that
%      is numerically rank deficient makes the DOC minimiser non-unique, which
%      makes the fit meaningless: it is being scored against an arbitrary
%      member of a solution set.
%
%      This is not hypothetical.  Both runs start from Q0 = I (cond = 1), so
%      the optimiser WALKS to the bad region.  Controlled comparison on
%      subject 4 / speed 1 / leg 1 (11 samples x 2 trials), same data, same
%      start, only cond_max differing:
%
%          cond_max = Inf,  40 iters:  E = 8.49   cond(Q) = 9.6e17
%          cond_max = 1e4, 150 iters:  E = 8.44   cond(Q) = 188
%
%      Identical fit quality, fifteen orders of magnitude of conditioning.
%      The ill-conditioning buys nothing, so there is no accuracy argument
%      for leaving it unconstrained.
%
%   Both are fixed by the pair (shift, trace constraint):
%
%       lambda_min(Q) >= eps_shift              (from the +eps_shift*I)
%       lambda_max(Q) <= trace(Q) = n           (all eigenvalues positive)
%       =>  cond(Q) <= n / eps_shift = cond_max     -- an exact guarantee
%
%   The trace constraint is imposed by the caller, not here; see
%   QP_TRACE_CONSTRAINT and QP_IO_FMINCON_SEARCH.  It must be an EQUALITY.
%   With trace(Q) <= n the entire gauge ray below the boundary stays feasible
%   at identical objective value, so the flat direction survives inside the
%   feasible set and nothing pushes the iterate to the boundary.  Only the
%   equality collapses the ray to a point.
%
%   Note the guarantee needs both halves.  The shift alone bounds lambda_min
%   but leaves lambda_max free; the trace constraint alone bounds lambda_max
%   but leaves lambda_min free.
%
%
%   ON THE GRADIENT
%   ---------------
%   KKT_SENSITIVITY needs no change.  It uses Q only in the cost Hessian
%   H = W'*Q*W, and L only in the dQ/dL term, and dQ/dL is unaffected by a
%   constant shift.  So the shifted Q and the UNSHIFTED L must be passed
%   together -- which is what this function returns.
%
%   Inputs:
%     theta    - (n*(n+1)/2 + n x 1) parameter vector [q; l]
%     n        - number of muscles
%     cond_max - (optional) condition-number ceiling, default 1e4; Inf disables
%
%   Outputs:
%     Q         - (n x n) symmetric, eigenvalues in [eps_shift, ...]
%     l         - (n x 1) linear weight vector
%     L         - (n x n) lower-triangular factor, WITHOUT the shift
%     eps_shift - the applied shift, n/cond_max (0 when cond_max is Inf)

if nargin < 3 || isempty(cond_max)
    cond_max = 1e4;
end

% Checked in two steps: a non-scalar cond_max cannot be interpolated with %g,
% and attempting it makes error() itself throw MATLAB:error:nonScalarInput,
% masking the identifier callers match on.
if ~isscalar(cond_max) || ~isreal(cond_max)
    error('theta_to_Ql:badCondMax', ...
          'cond_max must be a real scalar (or Inf); got a %s.', ...
          mat2str(size(cond_max)));
end

if cond_max <= 0 || isnan(cond_max)
    error('theta_to_Ql:badCondMax', ...
          'cond_max must be positive, got %g.', cond_max);
end

nq = n*(n+1)/2;

if numel(theta) ~= nq + n
    error('theta_to_Ql:badTheta', ...
          'theta has %d entries; expected n*(n+1)/2 + n = %d for n = %d.', ...
          numel(theta), nq + n, n);
end

if isinf(cond_max)
    eps_shift = 0;
else
    eps_shift = n / cond_max;
end

[Q, L] = chol_vec_to_Q(theta(1:nq), n);
Q = Q + eps_shift * eye(n);
Q = (Q + Q') / 2;                  % kill round-off asymmetry
l = reshape(theta(nq+1:end), n, 1);

end
