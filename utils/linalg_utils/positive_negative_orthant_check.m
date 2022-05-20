function [flag, varargout] = positive_negative_orthant_check(V)
%POSITIVE_NEGATIVE_ORTHANT_CHECK checks whether a vector is in the positive
%or negative orthant.
%
%   flag = POSITIVE_NEGATIVE_ORTHANT_CHECK(V) checks V.
%   [flag, Vp] = POSITIVE_NEGATIVE_ORTHANT_CHECK(V) checks V and returns
%   its positive orthant representation if yes.

% Check if all elements are positive of negative
flag = all(V >= 0) || all(V < 0);

% If additional argument are required
if nargout > 1
    % Positive orthant representation
    if flag
        % If it can be represent
        if all(V >= 0)
            Vp = V;
        else
            Vp = -V;
        end        
        varargout{1} = Vp;
    else
        varargout{1} = V;
    end
end
end