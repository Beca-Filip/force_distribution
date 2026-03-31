function [E, varargout] = QP_IO_inner_loop(theta, data, vars, model, sample_list, trial_list, speed_list, leg_list)
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

n = size(vars.variables.f, 1);
[Q, l, L] = theta_to_Ql(theta, n);

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
    % dE/dtheta = (1/(E*N)) * sum_{i,k,t} (Fout - Fref)[i,k,t] * dFout[i,p,k,t]
    %
    % Fout, Fref : [n x 1 x nsamples x ntrials x 1 x 1]
    % dFout      : [n x ntheta x nsamples x ntrials x 1 x 1]
    %
    % Reshape to [n x M] and [n x ntheta x M] where M = nsamples*ntrials,
    % then contract the n and M dimensions.

    ntheta = numel(theta);
    N      = numel(Fref);
    ns     = length(sample_list);
    nt     = length(trial_list);
    M      = ns * nt;

    r  = reshape(Fout(:,1,:,:,1,1) - Fref(:,1,:,:,1,1), n, M);  % [n x M]
    dF = reshape(dFout(:,:,:,:,1,1), n, ntheta, M);               % [n x ntheta x M]

    % dE[p] = (1/(E*N)) * sum_{i,m} r[i,m] * dF[i,p,m]
    %       = (1/(E*N)) * (dF_mat * r_vec)
    % where dF_mat = [ntheta x n*M], r_vec = [n*M x 1] (same column-major order)
    dF_mat = reshape(permute(dF, [2, 1, 3]), ntheta, n*M);  % [ntheta x n*M]
    r_vec  = r(:);                                            % [n*M x 1]

    varargout{1} = dF_mat * r_vec / (E * N);   % [ntheta x 1]
end

end

% -------------------------------------------------------------------------

function [Q, l, L] = theta_to_Ql(theta, n)
nq = n*(n+1)/2;
[Q, L] = chol_vec_to_Q(theta(1:nq), n);
l = theta(nq+1:end);
end
