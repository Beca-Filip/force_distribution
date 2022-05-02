function [Es, alphas] = IO_sort_rmses(E, alpha)
% Sort
[Es, Inds] = sort(E);
alphas = alpha(Inds, :);
end