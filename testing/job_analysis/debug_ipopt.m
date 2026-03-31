% Debugging CASADI ipopt


fprintf("-----------------------------------------------------------------\n");
% Check force values
debug_fmin = model.debug.value(vars.parameters.fmin);
debug_f = model.debug.value(vars.variables.f);
debug_fmax = model.debug.value(vars.parameters.fmax);
fprintf("Search variable lower bound values:\n");
fprintf("%.2f, ", debug_fmin);
fprintf("\n");
fprintf("Search variable values:\n");
fprintf("%.2f, ", debug_f);
fprintf("\n");
fprintf("Search variable upper bound values:\n");
fprintf("%.2f, ", debug_fmax);
fprintf("\n");

fprintf("-----------------------------------------------------------------\n");
% Check function values
debug_Jset = model.debug.value(vars.functions.Jset);
fprintf("Objective function values:\n");
fprintf("%.2e, ", debug_Jset);
fprintf("\n");

% Check function gradient values
debug_dJset = model.debug.value(jacobian(vars.functions.Jset, vars.variables.f)');

fprintf("Objective function gradient norms:\n");
fprintf("%.2e, ", sqrt(sum(debug_dJset.^2, 1)));
fprintf("\n");

% Check function hessian values
debug_ddJset = zeros(size(vars.variables.f, 1), size(vars.variables.f, 1), length(vars.functions.Jset));
debug_detddJset = zeros(1, length(vars.functions.Jset));
for ff = 1 : length(vars.functions.Jset)
    debug_ddJset(:, :, ff) = model.debug.value(hessian(vars.functions.Jset(ff), vars.variables.f));
    debug_detddJset(ff) = det(debug_ddJset(:, :, ff));
end

fprintf("Objective function hessian determinants:\n");
fprintf("%.2e, ", debug_detddJset);
fprintf("\n");

fprintf("-----------------------------------------------------------------\n");
% Check inequality constraint values
debug_c = model.debug.value(vars.functions.c);
fprintf("Are there any violated inequality constraints?: %d \n", any(debug_c > 0));
fprintf("Inequality constraint infeasibility (i.e. maximum constraint violation): %.2e \n", max(debug_c));
fprintf("Inequality constraint infeasible indices: \n")
fprintf("%d, ", find(debug_c > 0));
fprintf("\n");

% Check inequality constraint bounds
debug_fmin = model.debug.value(vars.parameters.fmin);
debug_fmax = model.debug.value(vars.parameters.fmax);
fprintf("Inequality constraint bound gap: \n");
fprintf("%.2e, ", debug_fmax - debug_fmin);
fprintf("\n");

fprintf("-----------------------------------------------------------------\n");
% Check equality constraint values
epsilon = 1e-3;
debug_ceq = model.debug.value(vars.functions.ceq);
fprintf("Equality constraint infeasibility: %.2e \n", max(abs(debug_ceq)));
fprintf("Equality constraint infeasible indices: \n")
fprintf("%d, ", find(abs(debug_ceq) > epsilon));
fprintf("\n");

