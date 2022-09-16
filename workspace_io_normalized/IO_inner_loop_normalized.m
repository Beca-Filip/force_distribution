function E = IO_inner_loop_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list)

% Initialize error function
E = zeros(size(data.f, 1), length(trial_list), length(speed_list), length(leg_list));
% Initialize total number of elemenst
totnum = 0;

% Get data forces
Fref = data.f(:,:,sample_list, trial_list, speed_list, leg_list);
% Get prediction forces
Fout = DO_subroutine_normalized(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list);

% Mean RMSE
% For each demonstration
for itr = 1 : length(trial_list)
    for isp = 1 : length(speed_list)
        for ile = 1 : length(leg_list)
            % For each muscle
            for mm = 1 : size(data.f, 1)                
                E(mm, itr, isp, ile) = rmse(data.f(mm, :, :, trial_list(itr), speed_list(isp), leg_list(ile)), Fout(mm, :, :, itr, isp, ile));
            end
        end
    end
end
E = mean(E(:));
end