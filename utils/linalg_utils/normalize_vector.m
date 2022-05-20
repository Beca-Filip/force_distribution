function Vn = normalize_vector(V, varargin)
%NORMALIZE_VECTOR normalizes vector to have norm 1, given a norm.
%
%   Vn = NORMALIZE_VECTOR(V) normalizes with norm 2.
%   Vn = NORMALIZE_VECTOR(V, n) normalizes with norm n.

% If additional argument passed
if nargin > 1
    n = varargin{1};
else
    n = 2;
end

% Check if n is infinity, return infinity norm
if n == inf
    Vn = V ./ max(abs(V));
    return
end

% Check if n is minus infinity, return minus infinity norm
if n == -inf
    Vn = V ./ min(abs(V));
    return
end

% If n is less than one or isn't equal to itself
if n < 1 || (n~=round(n))
    error('n should be a positive integer >= 1, inf or -inf.');
end

% For other cases
Vn = V ./ norm(V, n);
return
end