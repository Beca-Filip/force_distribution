function E = debug_IO_inner_loop_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list)

% Initialize error function
E = 0;
% Initialize total number of elemenst
totnum = 0;

% Get data forces
Fref = data.f(:,:,sample_list, trial_list, speed_list, leg_list);
% Get prediction forces
Fout = debug_DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% RMSE
E = rmse(Fref, Fout);

end