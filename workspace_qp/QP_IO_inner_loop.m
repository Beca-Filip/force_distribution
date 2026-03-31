function E = QP_IO_inner_loop(theta, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%QP_IO_INNER_LOOP Compute RMSE between QP-predicted and reference forces.
%
%   E = QP_IO_INNER_LOOP(theta, data, vars, model, sample_list, trial_list,
%                        speed_list, leg_list)
%
%   theta is a vector of length n*(n+1)/2 + n encoding the QP cost weights:
%     theta(1 : n*(n+1)/2)     -> lower-triangular entries of L, where Q = L*L'
%     theta(n*(n+1)/2+1 : end) -> l (n×1 linear weight vector)
%
%   The cost is: 0.5 * f_norm' * Q * f_norm + l' * f_norm
%   where f_norm = F_invcov * (f - f_mean).

n = size(vars.variables.f, 1);
[Q, l] = theta_to_Ql(theta, n);

% Get reference forces
Fref = data.f(:, :, sample_list, trial_list, speed_list, leg_list);

% Get QP-predicted forces
Fout = QP_subroutine(Q, l, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% RMSE
E = rmse(Fref, Fout);

end

% -------------------------------------------------------------------------

function [Q, l] = theta_to_Ql(theta, n)
nq = n*(n+1)/2;
[Q, ~] = chol_vec_to_Q(theta(1:nq), n);
l = theta(nq+1:end);
end
