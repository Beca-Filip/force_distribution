function syn_info = build_synergy_info(muscle_names, synergy_groups)
%BUILD_SYNERGY_INFO  Map named synergy groups to muscle indices and build
%                    the reduced block-diagonal Cholesky parameterisation.
%
%   syn_info = BUILD_SYNERGY_INFO(muscle_names, synergy_groups)
%
%   For each synergy group k of size m_k, the corresponding block L_k is a
%   m_k×m_k lower-triangular matrix with m_k*(m_k+1)/2 free parameters.
%   Q = blkdiag(L_1*L_1', ..., L_K*L_K', diag(l_s.^2))  is PSD by construction
%   and satisfies Q_{ij} = 0 whenever muscles i and j belong to different groups.
%   Muscles not assigned to any group are treated as singletons (one diagonal
%   free parameter each).
%
%   Each muscle is assigned to at most one group.  When multiple patterns in
%   the same group match a muscle that is a common issue, first-match-wins
%   within a group; across groups, the earlier group wins.
%
%   Inputs:
%     muscle_names   - {n x 1} cell array of muscle name strings
%     synergy_groups - {K x 1} cell array; each entry is a {1 x ?} cell of
%                      name patterns (partial, case-insensitive)
%
%   Output:  syn_info struct with fields
%     .n                 - total number of muscles
%     .n_groups          - K
%     .group_indices     - {K x 1} sorted muscle indices for each group
%     .singleton_indices - row vector of unassigned muscle indices
%     .n_full_q          - n*(n+1)/2  (length of the full q vector)
%     .n_theta_q         - number of free synergy q parameters
%     .n_theta           - n_theta_q + n  (full synergy theta length)
%     .full_q_indices    - [n_theta_q x 1]  position of each syn q param in
%                          the full q vector (1-indexed, column-major tril order)
%     .diag_mask         - logical [n_theta_q x 1], true for diagonal L entries
%                          (these must be >= 0 in fmincon lower bounds)

n  = numel(muscle_names);
K  = numel(synergy_groups);

% ---- Assign muscles to groups (first-group-first-match wins) -----------
assignment = zeros(n, 1);
for k = 1:K
    patterns = synergy_groups{k};
    for pi = 1:numel(patterns)
        pat = lower(patterns{pi});
        for mi = 1:n
            if assignment(mi) == 0 && contains(lower(muscle_names{mi}), pat)
                assignment(mi) = k;
            end
        end
    end
end

group_indices     = cell(K, 1);
for k = 1:K
    group_indices{k} = find(assignment == k)';   % row vector, sorted ascending
end
singleton_indices = find(assignment == 0)';

% ---- Index look-up table: (row, col) -> position in full q vector ------
% chol_vec_to_Q uses L(tril(true(n))) = q, which fills in column-major order.
[full_r, full_c] = find(tril(true(n)));
n_full_q         = numel(full_r);

% Build reverse map: q_pos_map(r, c) = position in full q vector
q_pos_map = zeros(n, n);
for k = 1:n_full_q
    q_pos_map(full_r(k), full_c(k)) = k;
end

% ---- Collect synergy q parameter indices (groups, then singletons) -----
full_q_idx_list = [];
diag_mask_list  = [];

for k = 1:K
    gi = group_indices{k};
    if isempty(gi), continue; end
    % Lower-triangular pairs within the group (column-major within-block order)
    for col_local = 1:numel(gi)
        for row_local = col_local:numel(gi)   % row >= col in lower-tri
            row_g = gi(row_local);
            col_g = gi(col_local);
            full_q_idx_list(end+1) = q_pos_map(row_g, col_g);  %#ok<AGROW>
            diag_mask_list(end+1)  = (row_g == col_g);          %#ok<AGROW>
        end
    end
end

for si = singleton_indices
    full_q_idx_list(end+1) = q_pos_map(si, si);  %#ok<AGROW>
    diag_mask_list(end+1)  = true;                %#ok<AGROW>
end

full_q_indices = full_q_idx_list(:);
diag_mask      = logical(diag_mask_list(:));
n_theta_q      = numel(full_q_indices);

% ---- Assemble output struct ---------------------------------------------
syn_info.n                 = n;
syn_info.n_groups          = K;
syn_info.group_indices     = group_indices;
syn_info.singleton_indices = singleton_indices;
syn_info.n_full_q          = n_full_q;
syn_info.n_theta_q         = n_theta_q;
syn_info.n_theta           = n_theta_q + n;
syn_info.full_q_indices    = full_q_indices;
syn_info.diag_mask         = diag_mask;

end
