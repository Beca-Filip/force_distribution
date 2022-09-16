function pcc = pearson_correlation_coefficient(x, y)
% Calculates the pearson correlation coefficient for a population of scalar
% random variables.

mu_x = mean(x);
mu_y = mean(y);

sig_x = std(x);
sig_y = std(y);

pcc = mean( (x-mu_x) .* (y-mu_y ) ) ./ (sig_x * sig_y);
end