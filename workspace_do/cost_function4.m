function J = cost_function4(f)
% APPROXIMATION OF MAXIMUM (INFINITY NORM)

% SOFTMAX
% % define scaling factor to avoid overload
% s = 100;
% 
% n = length(f);
% J = s*log(sum(exp(f ./ s)));

% 8-NORM
n = length(f);
J = (sum(f.^8) / n).^(1/8);

end