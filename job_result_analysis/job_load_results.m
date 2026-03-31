function results = job_load_results(resdir)
%JOB_LOAD_RESULTS loads the results of IOC local search. Assumes the
%directory string ends with slash or backslash.

% Append backslash
if ~(resdir.endsWith("/") || resdir.endsWith("\"))
    resdir = resdir.append("\");
end

% List of leg, speed, and phase
ls_l = 1 : 2;
ls_s = 1 : 2;
ls_p = 1 : 2;

warning("off");

% Output parameters
results = struct();
% Counter
counterVar = 1;

for l = ls_l
    for s = ls_s
        for p = ls_p
            % Load the corresponding file
            load(sprintf("%sl%ds%dp%d.mat", resdir, l, s, p));
            
            
            % Check if there exists a trial where simultaneously:
            %   - difference between any local alpha is greater than 1e-3
            %   - difference between any local error is greater than 5e-1
            %   newton
            if any( max( abs( diff(alpha_loc, 1, 1) ) , [], 1) > 1e-3) && any( max( abs( diff(E_loc, 1, 1) ) , [], 1)  > 5e-1)
                fprintf("l%ds%dp%d\n", l, s, p);
                alpha_loc
            end
            
            % Sort the local errors
            [~, errorSort] = sort(E_loc, 'ascend');
                        
            % Get the parametrization with least error
            results(counterVar).alpha = alpha_loc(errorSort(1), :);
            results(counterVar).error = E_loc(errorSort(1));
            
            % Store the leg, speed, trial list, and sample_list
            results(counterVar).leg = leg_list;
            results(counterVar).speed = speed_list;
            results(counterVar).trial_list = trial_list;
            results(counterVar).sample_list = sample_list;
            
            % Counter var increment
            counterVar = counterVar + 1;
        end
    end
end
warning("on");
end