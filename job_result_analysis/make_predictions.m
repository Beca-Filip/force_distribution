function data_ioc = make_predictions(data, results, model, vars)
%MAKE_PREDICTIONS creates a new f_pred field in the data structure data.

% Copy the data structure
data_ioc = data;
data_ioc.f_pred = zeros(size(data.f));

% For each result element do the DOC
for ii = 1 : length(results)
    Fout = DO_subroutine_normalized(results(ii).alpha, data, vars, model, ...
           results(ii).sample_list, results(ii).trial_list, results(ii).speed, results(ii).leg);
    
    data_ioc.f_pred(:, :, results(ii).sample_list, results(ii).trial_list, results(ii).speed, results(ii).leg) = Fout;


    % Make predictions for each cost function sparately
    for cf = 1 : length(vars.parameters.alpha)

        % Set cost function parameters
        alpha = zeros(1, length(vars.parameters.alpha));
        alpha(cf) = 1;

        % Make predictions
        Fout = DO_subroutine_normalized(alpha, data, vars, model, ...
               results(ii).sample_list, results(ii).trial_list, results(ii).speed, results(ii).leg);

        % Make predictions for each CF
        data_ioc.cf(cf).f_pred(:, :, results(ii).sample_list, results(ii).trial_list, results(ii).speed, results(ii).leg) = Fout;

    end
end

end