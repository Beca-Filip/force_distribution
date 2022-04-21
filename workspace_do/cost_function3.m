function J = cost_function3(f)

n = length(f);
J = (sum(f.^3) / n).^(1/3);

end