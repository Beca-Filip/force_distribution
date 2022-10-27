function [err, alpha] = debug_IO_grid_search_normalized(Ngrid, data, vars, model, sample_list, trial_list, speed_list, leg_list,debug_indices)
%IO_GRID_SEARCH_NORMALIZED performs a grid search over the normalized cost function
%parametrization.
%
%   [alpha, err] = IO_GRID_SEARCH_NORMALIZED(Ngrid, data, vars, model, sample_list, trial_list, speed_list, leg_list)
%   

% Cost function parametrization dimension and affine dimension
d_simplex = length(vars.parameters.J_max);
aff_d_simplex = d_simplex - 1;


% Determine the number of grid points to be used
m = determine_simplex_grid_partition(aff_d_simplex, Ngrid);

% Get the probability simplex gridpoints
alpha = prob_simplex_ndim(aff_d_simplex, m);

% Get the error vector
err = zeros(size(alpha, 1), 1);

% For each parametrization in the grid
for ii = debug_indices
    fprintf(strcat("ii=%d\n", nz_vector_elements_string(alpha(ii, :))), ii);
    % Calculate current rmse
    err(ii) = debug_IO_inner_loop_normalized(alpha(ii, :).', data, vars, model, sample_list, trial_list, speed_list, leg_list);

end

end