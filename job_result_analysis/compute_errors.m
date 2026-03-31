function results = compute_errors(data, results)
%COMPUTE_ERRORS computes all the errors and similarity measures that we
%want.
%Notably, RMSE per muscle and per trial, means of RMSE across trials per
%muscle, std of RMSE across trials per muscle, mean of (mean of RMSE across
%trials) across muscles, std of (mean of RMSE across trials) across
%muscles.
%Also, total RMSE across muscles and trials.
%
%Notably, CC per muscle and per trial, means of CC across trials per
%muscle, std of CC across trials per muscle, mean of (mean of CC across
%trials) across muscles, std of (mean of CC across trials) across
%muscles.
%
%Also, total CC across muscles and trials.

% For all results instances
for ii = 1 : length(results)
    % Take out the reference forces
    fref = data.f(:, :, results(ii).sample_list, results(ii).trial_list, results(ii).speed, results(ii).leg);
    % Take out the predicted forces
    fpred = data.f_pred(:, :, results(ii).sample_list, results(ii).trial_list, results(ii).speed, results(ii).leg);
    
    % Compute RMSE across muscle trajectories and trials
    results(ii).musctraj_rmse.values = muscle_traj_and_trial_rmse(fref, fpred);
    
    % Compute mean and std across trials, and across trajectories
    results(ii).musctraj_rmse.mean_across_trials = mean(results(ii).musctraj_rmse.values, 2);
    results(ii).musctraj_rmse.std_across_trials = std(results(ii).musctraj_rmse.values, [], 2);
    
    % Compute mean and std across 
    results(ii).musctraj_rmse.mean_across_trials_across_muscles = ...
        mean(results(ii).musctraj_rmse.mean_across_trials);
    results(ii).musctraj_rmse.std_across_trials_across_muscles = ...
        std(results(ii).musctraj_rmse.std_across_trials, [], 1);
    
    % Compute total RMSE for whole traj
    results(ii).total_rmse = rmse(fref, fpred);
    
    
    % Compute CC across muscle trajectories and trials
    results(ii).musctraj_cc.values = muscle_traj_and_trial_cc(fref, fpred);
    
    % Compute mean and std across trials, and across trajectories
    results(ii).musctraj_cc.mean_across_trials = mean(results(ii).musctraj_cc.values, 2);
    results(ii).musctraj_cc.std_across_trials = std(results(ii).musctraj_cc.values, [], 2);
    
    % Compute mean and std across 
    results(ii).musctraj_cc.mean_across_trials_across_muscles = ...
        mean(results(ii).musctraj_cc.mean_across_trials);
    results(ii).musctraj_cc.std_across_trials_across_muscles = ...
        std(results(ii).musctraj_cc.std_across_trials, [], 1);
    
    % Compute total CC for whole traj
    results(ii).total_cc = corr2(fref(:), fpred(:));
    
    % For each cost function separately compute prediction errors
    for cf = 1 : length(data.cf)
        % Take out the predicted forces
        fpred = data.cf(cf).f_pred(:, :, results(ii).sample_list, results(ii).trial_list, results(ii).speed, results(ii).leg);
    
        
        % Compute RMSE across muscle trajectories and trials
        results(ii).cf(cf).musctraj_rmse.values = muscle_traj_and_trial_rmse(fref, fpred);

        % Compute mean and std across trials, and across trajectories
        results(ii).cf(cf).musctraj_rmse.mean_across_trials = mean(results(ii).cf(cf).musctraj_rmse.values, 2);
        results(ii).cf(cf).musctraj_rmse.std_across_trials = std(results(ii).cf(cf).musctraj_rmse.values, [], 2);

        % Compute mean and std across 
        results(ii).cf(cf).musctraj_rmse.mean_across_trials_across_muscles = ...
            mean(results(ii).cf(cf).musctraj_rmse.mean_across_trials);
        results(ii).cf(cf).musctraj_rmse.std_across_trials_across_muscles = ...
            std(results(ii).cf(cf).musctraj_rmse.std_across_trials, [], 1);

        % Compute total RMSE for whole traj
        results(ii).cf(cf).total_rmse = rmse(fref, fpred);


        % Compute CC across muscle trajectories and trials
        results(ii).cf(cf).musctraj_cc.values = muscle_traj_and_trial_cc(fref, fpred);

        % Compute mean and std across trials, and across trajectories
        results(ii).cf(cf).musctraj_cc.mean_across_trials = mean(results(ii).cf(cf).musctraj_cc.values, 2);
        results(ii).cf(cf).musctraj_cc.std_across_trials = std(results(ii).cf(cf).musctraj_cc.values, [], 2);

        % Compute mean and std across 
        results(ii).cf(cf).musctraj_cc.mean_across_trials_across_muscles = ...
            mean(results(ii).cf(cf).musctraj_cc.mean_across_trials);
        results(ii).cf(cf).musctraj_cc.std_across_trials_across_muscles = ...
            std(results(ii).cf(cf).musctraj_cc.std_across_trials, [], 1);

        % Compute total CC for whole traj
        results(ii).cf(cf).total_cc = corr2(fref(:), fpred(:));
    end
end



end