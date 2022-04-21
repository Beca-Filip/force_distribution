function J = cost_function16(f, f0, r)

maxapprox = @(x) log(sum(exp( x )));

fnormalisedf0 = f./f0;

n = length(f);
J = maxapprox( (exp(3.48 + 0.169*r) .* (100*fnormalisedf0).^(-0.5 - 0.036*r)).^(-1) );

end