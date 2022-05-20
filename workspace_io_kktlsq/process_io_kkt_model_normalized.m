function modelio = process_io_kkt_model_normalized(modelio, aithresh)
%PROCESS_IO_KKT_MODEL_NORMALIZED forms the Keshavaraz et al. (2011) inverse
%optimization formulation model using the data, stores it inside a struct
%which contains numerical matrices needed for an LSQ resolution.
%
%   modelio = PROCESS_IO_KKT_MODEL_NORMALIZED(modelio, aithresh)
%   aithresh ~ active inequality treshold

% Extract lists
sample_list = modelio.sample_list;
trial_list = modelio.trial_list;
speed_list = modelio.speed_list;
leg_list = modelio.leg_list;

% Prealocate active inequality index cell array
modelio.aic = cell(numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));
% Prealocate number of active inequalities
modelio.num_aic = zeros(numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));
% Prealocate active inequality gradients cell array
modelio.gradients_aic = cell(numel(sample_list), numel(trial_list), numel(speed_list), numel(leg_list));

% Count number of active inequality constraints inside loop
num_tot_aic = 0;
% Counters for the loop
cntsamples = 1;
cnttrials = 1;
cntspeeds = 1;
cntlegs = 1;
% Loop over trials
for trial = trial_list
    % Reset speed counter
    cntspeeds = 1;
    % Loop over speeds
    for speed = speed_list
        % Reset leg counter
        cntlegs = 1;
        % Loop over legs
        for leg = leg_list
            % Reset sample counter
            cntsamples = 1;
            % Loop over samples
            for k = sample_list
                
                % Determine active inequality constraints and store them
                % (active are those that are superior to threshold)
                aic = find(modelio.values_ic(:, cntsamples, cnttrials, cntspeeds, cntlegs) >= aithresh);
                modelio.aic{cntsamples, cnttrials, cntspeeds, cntlegs} = aic;

                % Get number of active inequality constraints
                naic = length(aic);
                modelio.num_aic(cntsamples, cnttrials, cntspeeds, cntlegs) = naic;
                % Update total number of active inequality constraints
                num_tot_aic = num_tot_aic + naic;
                
                % Extract active inequality gradients
                dgaic = modelio.gradients_ic(:, aic, cntsamples, cnttrials, cntspeeds, cntlegs);
                modelio.gradients_aic{cntsamples, cnttrials, cntspeeds, cntlegs} = dgaic;

                % Augment sample counter
                cntsamples = cntsamples + 1;
            end
            % Augment leg count
            cntlegs = cntlegs + 1;
        end
        % Augment speeds counter
        cntspeeds = cntspeeds + 1;
    end
    % Augment trial counter
    cnttrials = cnttrials + 1;
end


%%% Create the regressor matrix

% INDEXING OF THE REGRESSOR:
% Get the cumulative sum of the number of active inequality constraints
modelio.cum_num_aic = vertcat(0, cumsum(modelio.num_aic(:)));
% Calculate the total number of rows: number of data points *
% dimensionality of DO variable
num_rows = modelio.do_n * modelio.num_total_data_pts;
% Calculate the total number of columns: number of cost function parameters
% + number of equality constraint multipliers + number of active inequality
% constraint multipliers
num_cf_parameters = modelio.num_cf;
num_ec_multipliers = modelio.num_ec * modelio.num_total_data_pts;
num_aic_multipliers = num_tot_aic;
num_cols = num_cf_parameters + num_ec_multipliers + num_aic_multipliers;


% Store total number of active inequality constraints
modelio.num_tot_ec = num_ec_multipliers;
modelio.num_tot_aic = num_tot_aic;
% Store the numbers of rows and columns
modelio.num_rows_regressor = num_rows;
modelio.num_cols_regressor = num_cols;

% Prealocate total regressor matrix
modelio.regressor = zeros(num_rows, num_cols);
% Form block matrix
% Counters for the loop
cntsamples = 1;
cnttrials = 1;
cntspeeds = 1;
cntlegs = 1;
% Loop over trials
for trial = trial_list
    % Reset speed counter
    cntspeeds = 1;
    % Loop over speeds
    for speed = speed_list
        % Reset leg counter
        cntlegs = 1;
        % Loop over legs
        for leg = leg_list
            % Reset sample counter
            cntsamples = 1;
            % Loop over samples
            for k = sample_list
                
                % Get indices for current cost function parameters
                [ir, ic] = indexing_function(modelio, 'cf_parameters', k, trial, speed, leg);
                % Assign the value to the the matrix
                modelio.regressor(ir, ic) = modelio.gradients_cf(:, :, cntsamples, cnttrials, cntspeeds, cntlegs);
                
                % Get indices for current equality constraint multipliers
                [ir, ic] = indexing_function(modelio, 'ec_multipliers', k, trial, speed, leg);
                % Assign the value to the the matrix
                modelio.regressor(ir, ic) = modelio.gradients_ec(:, :, cntsamples, cnttrials, cntspeeds, cntlegs);
                
                % Get indices for current equality constraint multipliers
                [ir, ic] = indexing_function(modelio, 'aic_multipliers', k, trial, speed, leg);
                % Assign the value to the the matrix
                modelio.regressor(ir, ic) = modelio.gradients_aic{cntsamples, cnttrials, cntspeeds, cntlegs};
                                
                % Augment sample counter
                cntsamples = cntsamples + 1;
            end
            % Augment leg count
            cntlegs = cntlegs + 1;
        end
        % Augment speeds counter
        cntspeeds = cntspeeds + 1;
    end
    % Augment trial counter
    cnttrials = cnttrials + 1;
end

%%% CREATE CONSTRAINT MATRICES

% Prealocate total lower bound vector
modelio.lower_bound = cat(1, zeros(modelio.num_cf, 1), -inf*ones(modelio.num_tot_ec, 1), zeros(modelio.num_tot_aic, 1));

% Prealocate normalization vector for cf weights
modelio.normalization_cf_Aeq = zeros(1, modelio.num_cols_regressor);
modelio.normalization_cf_Aeq(1:modelio.num_cf) = 1;
modelio.normalization_cf_beq = 1;

end