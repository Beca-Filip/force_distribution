function E = IO_inner_loop_normalized_variable_phase(alpha, data, vars, model, sample_cell_list, trial_list, speed_list, leg_list)

% Number of trajectories to compute
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

% Initialize error function
E = zeros(ntrials, nspeeds, nlegs);

% Initialize total number of elemenst
totnum = 0;

% Get prediction forces
Fout = DO_subroutine_normalized_variable_phase(alpha, data, vars, model, sample_cell_list, trial_list, speed_list, leg_list);

% RMSE
for trial = 1 : ntrials
    for speed = 1 : nspeeds
        for leg = 1 : nlegs
            E(trial, speed, leg) = rmse(data.f(:, :, sample_cell_list{trial}, trial, speed, leg),...
                                        Fout{trial, speed, leg});
            
        end
    end
end

% RMSE    
E = mean(E);

end