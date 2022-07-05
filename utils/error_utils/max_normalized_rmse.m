function e = max_normalized_rmse(a,b)
%MAX_NORMALIZED_RMSE calculates the max-normalized mean square error of two
%vector signals of the same size.
%
%   e = MAX_NORMALIZED_RMSE(a,b)
%   
%   Inputs:
%   a, b ~ multi-dimensional array. the comparison occurs along the 1st
%   dimension.
%   

e = sqrt(sum( abs((a-b) ./ max(a, [], 2:ndims(a))), 'all') ./ numel(a));

end