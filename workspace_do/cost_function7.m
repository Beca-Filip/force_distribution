function J = cost_function7(f, fmin, fmax)

n = length(f);
J = (sum( ((f-fmin)./fmax).^3 ) / n).^(1/3);

end