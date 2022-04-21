function J = cost_function12(f, pcsa)
% APPROXIMATION OF MAXIMUM (INFINITY NORM)

% add scaling factor to avoid infinity
s = 1000;

n = length(f);
J = s*log(sum(exp( f./pcsa./s)));

end