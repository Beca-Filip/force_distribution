function [E, varargout] = QP_synergy_IO_inner_loop(theta_syn, syn_info, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%QP_SYNERGY_IO_INNER_LOOP  RMSE (and gradient) for synergy-constrained QP IO.
%
%   E        = QP_SYNERGY_IO_INNER_LOOP(theta_syn, syn_info, ...)
%   [E, dE]  = QP_SYNERGY_IO_INNER_LOOP(theta_syn, syn_info, ...)
%
%   Behaves identically to QP_IO_INNER_LOOP but operates in the reduced
%   synergy parameter space defined by syn_info (see build_synergy_info).
%
%   theta_syn has syn_info.n_theta components:
%     theta_syn(1 : n_theta_q)     -> free block-diagonal L entries
%     theta_syn(n_theta_q+1 : end) -> l (linear weight vector)
%
%   Gradient chain rule:
%     dE/d(theta_syn) = [ dE/d(q_full)(full_q_indices) ]
%                       [ dE/d(l)                       ]
%   where dE/d(theta_full) is the full gradient from QP_subroutine/kkt_sensitivity
%   and the q part is simply selected by index (the Jacobian d(q_full)/d(q_syn)
%   is a 0/1 selection matrix).

n         = syn_info.n;
n_theta_q = syn_info.n_theta_q;

[Q, L, q_full, l] = synergy_theta_to_Q(theta_syn, syn_info);

Fref      = data.f(:, :, sample_list, trial_list, speed_list, leg_list);
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
    % dFout : [n x n_theta_full x nsamples x ntrials x 1 x 1]
    % where n_theta_full = n_full_q + n
    n_theta_full = syn_info.n_full_q + n;
    N            = numel(Fref);
    ns           = length(sample_list);
    nt           = length(trial_list);
    M            = ns * nt;

    r      = reshape(Fout(:,1,:,:,1,1) - Fref(:,1,:,:,1,1), n, M);    % [n x M]
    dF     = reshape(dFout(:,:,:,:,1,1), n, n_theta_full, M);           % [n x n_theta_full x M]

    dF_mat    = reshape(permute(dF, [2, 1, 3]), n_theta_full, n*M);    % [n_theta_full x n*M]
    dE_full   = dF_mat * r(:) / (E * N);                               % [n_theta_full x 1]

    % Chain rule: select synergy q entries, keep l unchanged
    dE_q_syn  = dE_full(syn_info.full_q_indices);     % [n_theta_q x 1]
    dE_l      = dE_full(syn_info.n_full_q+1:end);     % [n x 1]

    varargout{1} = [dE_q_syn; dE_l];                  % [n_theta_syn x 1]
end

end
