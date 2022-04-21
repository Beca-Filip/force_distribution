function J = cost_function13(f, vmt)

n = length(f);
J = sqrt(sum((f.*vmt).^2) / n);

end