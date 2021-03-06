function [err, alpha] = IO_grid_search(Ngrid, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%IO_GRID_SEARCH performs a grid search over the cost function
%parametrization.
%
%   [alpha, err] = IO_GRID_SEARCH(Ngrid, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%   

% Cost function parametrization dimension and affine dimension
d_simplex = length(data.J_max);
aff_d_simplex = d_simplex - 1;


% Determine the number of grid points to be used
m = determine_simplex_grid_partition(aff_d_simplex, Ngrid);

% Get the probability simplex gridpoints
alpha = prob_simplex_ndim(aff_d_simplex, m);

% Get the error vector
err = zeros(size(alpha, 1), 1);

% For each parametrization in the grid
for ii = 1 : size(alpha, 1)
    
    % Calculate current rmse
    err(ii) = IO_inner_loop(alpha(ii, :).', data, vars, model, sample_list, trial_list, speed_list, leg_list);

end

end