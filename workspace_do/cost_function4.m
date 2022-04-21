function J = cost_function4(f)
% APPROXIMATION OF MAXIMUM (INFINITY NORM)

% define scaling factor to avoid overload
s = 100;

n = length(f);
J = s*log(sum(exp(f ./ s)));

end