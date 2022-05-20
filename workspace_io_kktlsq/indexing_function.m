function [indices_row, indices_col] = indexing_function(modelio, type_str, sample, trial, speed, leg)
%INDEXING_FUNCTION returns the row and column indices in the inverse 
%optimization regressor matrix of cost function parameters or
%equality/inequality multipliers of a particular data point, given the data
%point indices and the type of desired variable.
%   - cost function parameters when type_str is 'cf_parameters'
%   - equality constraint multipliers when type_str is 'ec_multipliers'
%   - active inequality constraint multipliers when type_str is
%   'aic_multipliers'


% Requested cost function parameters
if strcmp(type_str, 'cf_parameters')
    
    % If particular sample requested
    if ~isempty(sample) && ~isempty(trial) && ~isempty(speed) && ~isempty(leg)
        % Determine linear index
        lin_ind = data_point_linear_index(modelio, sample, trial, speed, leg);
        % The row index is determined as the number of DO variables times 
        % the linear index
        indices_row = modelio.do_n * (lin_ind-1) + 1 : modelio.do_n * lin_ind;
    % If not, give all samples
    else
        indices_row = 1 : modelio.num_rows_regressor;
    end
    % Columns are always the same: the first columns correspond to the cost
    % function parameters
    indices_col = 1 : modelio.num_cf;
    return
end

% Requested equality constraint multipliers
if strcmp(type_str, 'ec_multipliers')
        
    % If particular sample requested
    if ~isempty(sample) && ~isempty(trial) && ~isempty(speed) && ~isempty(leg)
        % From the linear index 
        lin_ind = data_point_linear_index(modelio, sample, trial, speed, leg);
        % Find the row indices as the number of DO variables times the 
        % linear index
        indices_row = modelio.do_n * (lin_ind-1) + 1 : modelio.do_n * lin_ind;
        
        % The column index is determined as the number of DO constraints times 
        % the linear index, plus the number of cost function parameters
        % which come before them
        indices_col = modelio.num_cf + (modelio.num_ec * (lin_ind-1) + 1 : modelio.num_ec * lin_ind);
    % If not, give all samples
    else
        indices_row = 1 : modelio.num_rows_regressor;
        indices_col = modelio.num_cf + 1 : modelio.num_cf + modelio.num_tot_ec;
    end
    return
end

% Requested inequality constraint multipliers
if strcmp(type_str, 'aic_multipliers')
        
    % If particular sample requested
    if ~isempty(sample) && ~isempty(trial) && ~isempty(speed) && ~isempty(leg)
        % From the linear index 
        lin_ind = data_point_linear_index(modelio, sample, trial, speed, leg);
        % Find the row indices as the number of DO variables times the 
        % linear index
        indices_row = modelio.do_n * (lin_ind-1) + 1 : modelio.do_n * lin_ind;
        
        % The column index is determined as the number of active DO 
        % constraints before this one, plus the number of cost function parameters
        % which come before them plus the total number of equality
        % constraint multipliers
        indices_col = modelio.num_cf + modelio.num_tot_ec + ...
                (modelio.cum_num_aic(lin_ind)+1:modelio.cum_num_aic(lin_ind+1));
    % If not, give all samples
    else
        indices_row = 1 : modelio.num_rows_regressor;
        indices_col = modelio.num_cf + modelio.num_tot_ec + 1 : modelio.num_cf + modelio.num_tot_ec + modelio.num_tot_aic;
    end    
    return
end

% If none of these output error
errmsg = "type_str should be one of the following three: cf_parameters, ec_multipliers, aic_multipliers.";
error(errmsg);

end