function nA = column_norm(A, varargin)
%COLUMN_NORM calculates the row-vector of column norms of a matrix.
%
%   nA = MULTIDIM_COLUMN_NORM(A)
%   nA = MULTIDIM_COLUMN_NORM(A, n)
%
%   Inputs:
%   A ~ matrix with size [n1, n2].
%   n ~ which norm to calculate (any positive integer or (+-) inf).
%   
%   Outputs:
%   nA ~ row-vector with size [1, n2] containing norms of columns of matrix
%   A.

% If additional argument passed
if nargin > 1
    n = varargin{1};
else
    n = 2;
end

% If n is not an integer greater or equal to one or isn't +- inf.
if ~( n == inf || n == -inf || (n >= 1 && n == round(n)) )
   error('n should be a positive integer >= 1, inf or -inf.');
end

% Get size
szA = size(A);

% Prealocate norm
nA = nan(1, szA(2));

% For each column
for cc = 1 : szA(2)
   
    % Calculate norm of the column
    nA(cc) = norm(A(:, cc), n);
    
end