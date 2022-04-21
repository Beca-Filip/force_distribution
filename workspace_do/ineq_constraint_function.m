function c = ineq_constraint_function(f, fmin, fmax)
c = [-f + fmin; f - fmax];
end