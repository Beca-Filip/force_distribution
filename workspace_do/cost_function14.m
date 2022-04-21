function J = cost_function14(f, M)

n = length(f);
J = sqrt(sum((f./M).^2) / n);

end