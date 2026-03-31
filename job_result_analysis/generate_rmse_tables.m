function generate_rmse_tables(results, dirname)
%GENERATE_RMSE_TABLES generates the RMSE and CC tables for a given subject 
%at a given speed. The table contains two large columns for the stance and 
%swing phase, each of which contains two subcolumns one for the RMSE and
%the second for the CC.
%
%   Example:
%   TABLE Leg -- Speed -- 
%           | RMSE PH1 [N] | CC PH1       |  RMSE PH2 [N] | CC PH2       |
% J_bilevel | --- +- --    | .-- +- .--   |  --- +- --    | .-- +- .--   |
% phi_1     | --- +- --    | .-- +- .--   |  --- +- --    | .-- +- .--   |
% [...]     | --- +- --    | .-- +- .--   |  --- +- --    | .-- +- .--   |
% phi_n     | --- +- --    | .-- +- .--   |  --- +- --    | .-- +- .--   |

% Define phase 1 and phase 2 sample lists
PH1 = 1:4:61;
PH2 = 61:4:101;

% Check dirname
if ~exist(dirname, 'dir')
    mkdir(dirname)
end

% Find speed leg pairs
speed_leg_pairs = [];
for rr = 1 : length(results)
    % current speed leg pair
    curr_sp_le_pair = [results(rr).speed, results(rr).leg];
    % if speed leg pairs is empty or
    % if current speed leg pair is not already a row in speed leg pairs
    if isempty(speed_leg_pairs) || ~any( ismember(speed_leg_pairs, curr_sp_le_pair, 'rows') )
        % add it
       speed_leg_pairs = [speed_leg_pairs; curr_sp_le_pair]; 
    end
end

% For each speed leg pair create a table
for ii = 1 : size(speed_leg_pairs, 1)
    
    % Table
    % Initialize to a table with 4 columns and rows corresponding to the
    % compound cost function and to each individual cost function
    numrows = 1 + length(results(1).cf);
    cfnames = arrayfun(@(ii) sprintf("$\\phi_{%d}$", ii), 1:length(results(1).cf), 'UniformOutput', true);
    rownames = horzcat("$J_{\rm bilevel}$", cfnames);
    tab = table('Size', [numrows, 4],...
        'VariableType', repmat("string", 1, 4),...
        'VariableNames', ["RMSE PH1 [N]", "CC PH1", "RMSE PH2 [N]", "CC PH2"], ...
        'RowNames', rownames);
    
    % Traverse each result instance
    for rr = 1 : length(results)
        
        % current speed leg pair
        curr_sp_le_pair = [results(rr).speed, results(rr).leg];
        % If result instance corresponds to given speed leg pair
        if isequal(speed_leg_pairs(ii, :), curr_sp_le_pair)
            
            % Check which phase
            PH = 0;
            % If sample list is similar to phase 1
            if isequal(results(rr).sample_list, PH1)
                PH = 1;
            else
                PH = 2;
            end
            
            % Create the string that describes the results of the bilevel
            rmse_string = sprintf("$%d \\pm %d$", round(results(rr).total_rmse),...
                                round(results(rr).musctraj_rmse.std_across_trials_across_muscles));
            cc_string = sprintf("$%.2f \\pm %.2f$", results(rr).total_cc,...
                                results(rr).musctraj_cc.std_across_trials_across_muscles);
            % CC make .xx instead of 0.xx
            cc_string = strrep(cc_string, "0.", ".");
            
            % Inscribe results for bilevel inside the table, choose
            % appropriate phase according to PH variable
            % Bilevel is 1st in table
            tab.(sprintf("RMSE PH%d [N]", PH))(1) = rmse_string;
            tab.(sprintf("CC PH%d", PH))(1) = cc_string;
            
            % Do the same for each cost function
            for cf = 1 : length(results(1).cf)
                
                % Create the string that describes the results of the said
                % cost function
                rmse_string = sprintf("$%d \\pm %d$", round(results(rr).cf(cf).total_rmse),...
                                round(results(rr).cf(cf).musctraj_rmse.std_across_trials_across_muscles));
                cc_string = sprintf("$%.2f \\pm %.2f$", results(rr).cf(cf).total_cc,...
                                results(rr).cf(cf).musctraj_cc.std_across_trials_across_muscles);
                % CC make .xx instead of 0.xx
                cc_string = strrep(cc_string, "0.", ".");
                % Inscribe results for teh said cost function inside the table, 
                % choose appropriate phase according to PH variable
                % The said cost function is (cf+1)-th in table
                tab.(sprintf("RMSE PH%d [N]", PH))(1+cf) = rmse_string;
                tab.(sprintf("CC PH%d", PH))(1+cf) = cc_string;
            end
            
        end
    end
    
    % Write the table
    filetype = ".xlsx";
    filename = sprintf("Speed_%d_leg_%d%s", speed_leg_pairs(ii, 1), ...
        speed_leg_pairs(ii, 2), filetype);
    filename = strcat(dirname, filename);
    
    writetable(tab, filename, 'WriteRowNames', true, 'FileType', 'Spreadsheet');
end
end