function arg1or2 = ternary_operator(cond, arg1, arg2)
%TERNARY_OPERATOR returns one of two inputs given a condition. Like the
%ternary operator in C.
%
%   arg1or2 = TERNARY_OPERATOR(cond, arg1, arg2) returns arg1 if cond is
%   true, and arg2 if cond is false.

if ~isscalar(cond) || isempty(cond)
    error('cond must be a non-empty scalar variable.');
end

if cond
    arg1or2 = arg1;
else
    arg1or2 = arg2;
end
end