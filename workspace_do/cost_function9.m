function J = cost_function9(f, pcsa)

n = length(f);
J = sum(f./pcsa) / n;

end