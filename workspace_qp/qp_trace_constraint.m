function [c, ceq, gc, gceq] = qp_trace_constraint(theta, n, cond_max)
%QP_TRACE_CONSTRAINT  Nonlinear equality fixing trace(Q) = n, for fmincon.
%
%   [c, ceq, gc, gceq] = QP_TRACE_CONSTRAINT(theta, n, cond_max)
%
%   Enforces trace(Q) = n where Q = L*L' + eps_shift*I is built by
%   THETA_TO_QL.  Since trace(L*L') = sum of squares of all entries of L,
%
%       trace(Q) = ||q||^2 + n*eps_shift  =  n
%       <=>  ceq(theta) = ||q||^2 - n*(1 - eps_shift) = 0
%
%   with eps_shift = n/cond_max.  There are no inequality constraints.
%
%   Together with the eps_shift*I floor on the eigenvalues this gives the
%   exact bound cond(Q) <= n/eps_shift = cond_max; see THETA_TO_QL for the
%   full argument and for why this must be an equality rather than "<=".
%
%   The constraint is a sphere of radius sqrt(n*(1-eps_shift)) in q, so the
%   feasible set is nonconvex -- but it is smooth, its gradient 2q never
%   vanishes on the sphere (the radius is positive), and fmincon's SQP and
%   interior-point algorithms handle it without trouble.  l is unconstrained,
%   which is what fixes the gauge: (Q,l) -> (cQ, cl) can no longer be
%   traversed, because c ~= 1 changes trace(Q).
%
%   Gradient convention: fmincon expects gceq to be [numel(theta) x n_ceq],
%   i.e. one COLUMN per equality constraint -- the transpose of the Jacobian.
%   Pass 'SpecifyConstraintGradient', true so these are used.
%
%   Inputs:
%     theta    - (n*(n+1)/2 + n x 1) parameter vector [q; l]
%     n        - number of muscles
%     cond_max - condition-number ceiling (Inf disables the shift)
%
%   Outputs:
%     c    - [] (no inequality constraints)
%     ceq  - scalar residual ||q||^2 - n*(1 - eps_shift)
%     gc   - [] (no inequality constraints)
%     gceq - (numel(theta) x 1) gradient of ceq

nq = n*(n+1)/2;

if isinf(cond_max)
    eps_shift = 0;
else
    eps_shift = n / cond_max;
end

q = theta(1:nq);

c   = [];
gc  = [];
ceq = q' * q - n * (1 - eps_shift);

gceq              = zeros(numel(theta), 1);
gceq(1:nq)        = 2 * q;

end
