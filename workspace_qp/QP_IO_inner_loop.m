function [E, varargout] = QP_IO_inner_loop(theta, data, vars, model, sample_list, trial_list, speed_list, leg_list, cond_max)
%QP_IO_INNER_LOOP  RMSE (and optionally its gradient) for the QP IO problem.
%
%   E        = QP_IO_INNER_LOOP(theta, ...)   objective only
%   [E, dE]  = QP_IO_INNER_LOOP(theta, ...)   objective + gradient
%
%   When called with two outputs the gradient dE/dtheta is computed
%   analytically via KKT sensitivity (see kkt_sensitivity.m) instead of
%   finite differences, enabling gradient-based IO search in fmincon.
%
%   theta is a vector of length n*(n+1)/2 + n:
%     theta(1 : n*(n+1)/2)     -> lower-triangular entries of L, Q = L*L'
%     theta(n*(n+1)/2+1 : end) -> l  (linear weight vector)
%
%   sample_list, trial_list, speed_list and leg_list may each be vectors.
%   Both E and dE are then pooled over every (sample, trial, speed, leg)
%   combination, so a single (Q, l) can be fitted jointly across speeds and
%   legs rather than one condition at a time.
%
%   cond_max is optional (default 1e4) and is forwarded to THETA_TO_QL, which
%   builds Q = L*L' + (n/cond_max)*I.  Pass Inf for the legacy unshifted
%   Q = L*L'.  The trace constraint that completes the conditioning guarantee
%   is imposed by QP_IO_FMINCON_SEARCH, not here; this function is just the
%   objective and stays well defined at any theta.

n = size(vars.variables.f, 1);

if nargin < 9 || isempty(cond_max)
    cond_max = 1e4;
end

% Q carries the eps_shift, L does not.  kkt_sensitivity needs exactly that
% pairing: Q enters the Hessian H = W'*Q*W, L enters dQ/dL, and dQ/dL is
% unaffected by a constant shift.
[Q, l, L] = theta_to_Ql(theta, n, cond_max);

Fref = data.f(:, :, sample_list, trial_list, speed_list, leg_list);

grad_flag = nargout > 1;

if grad_flag
    [Fout, dFout] = QP_subroutine(Q, l, data, vars, model, ...
        sample_list, trial_list, speed_list, leg_list, L);
else
    Fout = QP_subroutine(Q, l, data, vars, model, ...
        sample_list, trial_list, speed_list, leg_list);
end

E = rmse(Fref, Fout);

if grad_flag
    % dE/dtheta = (1/(E*N)) * sum_{i,m} (Fout - Fref)[i,m] * dFout[i,p,m]
    % where m ranges over every condition (sample, trial, speed, leg).
    %
    % Fout, Fref : [n x 1 x nsamples x ntrials x nspeeds x nlegs]
    % dFout      : [n x ntheta x nsamples x ntrials x nspeeds x nlegs]
    %
    % Reshape to [n x M] and [n x ntheta x M] where
    % M = nsamples*ntrials*nspeeds*nlegs, then contract the n and M dims.
    %
    % Dimension 2 of Fout/Fref is singleton, so collapsing dims 3..6 into a
    % single index m enumerates conditions in column-major order
    % (sample fastest, then trial, then speed, then leg) -- exactly the order
    % dFout's dims 3..6 collapse into.  The two are therefore aligned
    % element-for-element, and reshape itself errors out if a dimension is
    % ever mismatched.

    ntheta = numel(theta);
    N      = numel(Fref);
    ns     = length(sample_list);
    nt     = length(trial_list);
    nsp    = length(speed_list);
    nlg    = length(leg_list);
    M      = ns * nt * nsp * nlg;

    r  = reshape(Fout - Fref, n, M);           % [n x M]
    dF = reshape(dFout, n, ntheta, M);         % [n x ntheta x M]

    % dE[p] = (1/(E*N)) * sum_{i,m} r[i,m] * dF[i,p,m]
    %       = (1/(E*N)) * (dF_mat * r_vec)
    % where dF_mat = [ntheta x n*M], r_vec = [n*M x 1] (same column-major order)
    dF_mat = reshape(permute(dF, [2, 1, 3]), ntheta, n*M);  % [ntheta x n*M]
    r_vec  = r(:);                                            % [n*M x 1]

    varargout{1} = dF_mat * r_vec / (E * N);   % [ntheta x 1]
end

end
