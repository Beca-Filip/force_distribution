function Cab = vector_xcorr(a,b)
%VECTOR_XCORR calculates the vector cross correlation between two signals, 
%which are given in row form (meaning the signals are given in different 
%rows, with length equal to number of columns).
%
%   Cab = VECTOR_XCORR(a,b)

% Prealocate matrix of cross correlations
Cab = zeros(size(a, 1), 2 * size(a, 2) - 1);

% Calculate for each row
for ii = 1 : size(a, 1)
    Cab(ii, :) = xcorr(a(ii, :), b(ii, :));
end

end