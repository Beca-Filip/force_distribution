function [Q, L, q_full, l] = synergy_theta_to_Q(theta_syn, syn_info)
%SYNERGY_THETA_TO_Q  Expand reduced synergy parameters to Q = L*L'.
%
%   [Q, L, q_full, l] = SYNERGY_THETA_TO_Q(theta_syn, syn_info)
%
%   theta_syn is a vector of length syn_info.n_theta:
%     theta_syn(1 : n_theta_q)     -> free lower-triangular entries of the
%                                      block-diagonal L (one block per synergy
%                                      group plus one scalar per singleton)
%     theta_syn(n_theta_q+1 : end) -> l  (linear weight vector, n x 1)
%
%   The expansion is:
%     q_full = zeros(n_full_q, 1)
%     q_full(syn_info.full_q_indices) = theta_syn(1:n_theta_q)
%   which places the synergy parameters at their correct positions in the
%   full lower-triangular L vector and leaves all cross-synergy entries at 0.
%   Q = L*L' is then block-diagonal with the desired sparsity pattern.
%
%   Inputs:
%     theta_syn - [syn_info.n_theta x 1] synergy parameter vector
%     syn_info  - struct from build_synergy_info
%
%   Outputs:
%     Q      - [n x n] symmetric PSD weight matrix
%     L      - [n x n] lower-triangular Cholesky factor, Q = L*L'
%     q_full - [n_full_q x 1] full lower-tri L parameter vector (with zeros)
%     l      - [n x 1] linear weight vector

n_theta_q = syn_info.n_theta_q;

q_full                          = zeros(syn_info.n_full_q, 1);
q_full(syn_info.full_q_indices) = theta_syn(1:n_theta_q);
l                               = theta_syn(n_theta_q+1:end);

[Q, L] = chol_vec_to_Q(q_full, syn_info.n);

end
