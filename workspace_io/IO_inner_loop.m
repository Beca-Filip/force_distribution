function E = IO_inner_loop(alpha, data, vars, model, sample_list, trial_list, speed_list, leg_list)

% Initialize error function
E = 0;
% Initialize total number of elemenst
totnum = 0;

% Loop over trials
for trial = trial_list
    % Loop over speeds
    for speed = speed_list
        % Loop over legs
        for leg = leg_list
            % Loop over samples
            for k = sample_list

                % Set model parameters
                model = set_model_parameters(data, vars, model, k, trial, speed, leg);

                % Set model weights
                model = set_model_weights(alpha, vars, model);

                % Optimize
                sol = model.solve();

                % Take the solution forces
                f_opt = sol.value(vars.variables.f);

                % Calculate rmse
                E = (E + rmse(f_opt, data.f(:,:,k,trial,speed,leg)).^2 * numel(f_opt));
                totnum = totnum + numel(f_opt);
            end
        end
    end
end

E = sqrt(E / totnum);

end