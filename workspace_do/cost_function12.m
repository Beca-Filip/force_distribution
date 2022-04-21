function J = cost_function12(f, pcsa)
% APPROXIMATION OF MAXIMUM (INFINITY NORM)

% SOFTMAX
% % add scaling factor to avoid infinity
% s = 1000;
% 
% n = length(f);
% J = s*log(sum(exp( f./pcsa./s)));

% 8-NORM
n = length(f);
J = (sum( (f./pcsa).^8) / n).^(1/8);
end