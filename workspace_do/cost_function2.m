function J = cost_function2(f)

n = length(f);
J = sqrt(sum(f.^2) / n);

end