function J = cost_function5(f, fmin, fmax)

n = length(f);
J = sum( (f-fmin)./fmax ) / n;

end