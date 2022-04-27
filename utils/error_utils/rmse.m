function e = rmse(a,b)
%RMSE calculates the root mean square error of two signals of the same
%size.
%
%   e = RMSE(a,b)

e = sqrt(sum(sum((a-b).^2)) ./ numel(a));

end