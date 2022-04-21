function J = cost_function8(f, fmin, fmax)
% APPROXIMATION OF MAXIMUM (INFINITY NORM)

% SOFTMAX
n = length(f);
J = log(sum(exp( (f-fmin)./fmax )));

end