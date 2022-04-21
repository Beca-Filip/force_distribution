function J = cost_function10(f, pcsa)

n = length(f);
J = sqrt(sum((f./pcsa).^2) / n);

end