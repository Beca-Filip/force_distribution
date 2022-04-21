function J = cost_function11(f, pcsa)

n = length(f);
J = (sum((f./pcsa).^3) / n).^(1/3);

end