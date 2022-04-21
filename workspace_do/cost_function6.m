function J = cost_function6(f, fmin, fmax)

n = length(f);
J = sqrt(sum( ((f-fmin)./fmax).^2 ) / n);

end