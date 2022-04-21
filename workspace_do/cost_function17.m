function J = cost_function17(f, fmax)


n = length(f);
J = - sum(sqrt(1 - (f./fmax).^2)) / n;

end