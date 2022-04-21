function e = nrmse(a,b)
%RMSE calculates the normalized root mean square error of two signals of
%the same size.
%
%   e = RMSE(a,b)
%   a is the reference signal

e = sqrt(sum((a-b).^2 ./ a.^2, 'all') ./ numel(a));

end