% Exclude cost functions
cf_exclude = [16, 17];
% Create model
[model, vars] = form_casadi_model_normalized(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);

Fout = DO_subroutine_normalized(alpha_loc(3, :).',data,vars,model,sample_list,trial_list,speed_list,leg_list);

CC_loc = corr2(reshape(Fout, [], 1), reshape(data.f(:,:,sample_list,trial_list,speed_list,leg_list), [], 1));

%%
fprintf('E_loc = [%.0f, %.0f, %.0f, %.0f].\n', E_loc);
fprintf('CC_loc = [%.2f].\n', CC_loc);
%%
fprintf(strcat('alpha_loc = [', repmat('%.2f, ', 1, 14), '%.2f].\n'), alpha_loc(4, :));