function generate_alpha_tables(results, dirname)
%GENERATE_ALPHA_TABLES generates the cost function weight tables for a 
%given subject at a given speed, and with given leg.
%
%The table is grouped by leg, then by speed, then by phase.
%   Example:
%   TABLE Leg -- Speed -- 
%           | Leg 1                       |  Leg 2                       |
%           | Speed 1      | Speed 2      |  Speed 1      |  Speed 2     |
%           | P1   | P2    | P1   | P2    | P1    | P2    | P1    | P2   |  
% omega_1   | .xx  | .xx   | .xx  | .xx   | .xx   | .xx   | .xx   | .xx  |
% omega_2   | .xx  | .xx   | .xx  | .xx   | .xx   | .xx   | .xx   | .xx  |
% [...]     | .xx  | .xx   | .xx  | .xx   | .xx   | .xx   | .xx   | .xx  |
% omega_n   | .xx  | .xx   | .xx  | .xx   | .xx   | .xx   | .xx   | .xx  |

% Phase 1 sample list
PH1 = 1:4:61;

% Extract leg speed phase triads
leg_speed_phase = [];
% Extract corresponding bilevel gotten alpha
alpha_bilevel = [];
for rr = 1 : length(results)
    % Phase convert
    if isequal(PH1, results(rr).sample_list)
        PH = 1;
    else
        PH = 2;
    end
    
    % Current legspeedphase triad
    curr_lsp_triad = [results(rr).leg, results(rr).speed, PH];
    
    % Add current triad
    leg_speed_phase = [leg_speed_phase; curr_lsp_triad];
    
    % Add current alpha
    alpha_bilevel = horzcat(alpha_bilevel, reshape(results(rr).alpha, [], 1)); 
end

% Sort them (use trick to multiply columns by powers of ten, and to sort
% subsequently)
lsp_trick = leg_speed_phase(:, 1) * 1e2 + leg_speed_phase(:, 2) * 1e1 + leg_speed_phase(:, 3) * 1;
[~, lsp_row_order] = sort(lsp_trick, 'ascend');
% Reorder
leg_speed_phase = leg_speed_phase(lsp_row_order, :);


% Use the same sort indexing to reorder alpha columns
alpha_bilevel = alpha_bilevel(:, lsp_row_order);


%%% Table creation

% Number of rows is equal to the number of cost function parameters
numrows = size(alpha_bilevel, 1);
% The name of the rows corresponds to parameter names, which are omega
% indexed by the number of the cost function
rownames = arrayfun(@(ii) sprintf("$\\omega_{%d}$", ii), 1:numrows, 'UniformOutput', true);
% Number of columns is the number of results
numcols = length(results);


% Table constructor
tab = table('Size', [numrows, numcols],...
    'VariableType', repmat("string", 1, numcols),...
    'RowNames', rownames);

% Fill the table with the cost function parametrizations
% contained in alpha_bilevel
for rr = 1 : numrows
    for cc = 1 : numcols
        
        % String formating
        str_alpha = sprintf("%.2f", alpha_bilevel(rr, cc));
        % Replace the "0.xx" by ".xx" in the formatting 
        str_alpha = strrep(str_alpha, "0.", ".");
        
        tab{rr, cc} = str_alpha;
        
    end
end


% Write the table
filetype = ".xlsx";
filename = sprintf("recovered_cost_parameters%s", filetype);
filename = strcat(dirname, filename);

writetable(tab, filename, 'WriteRowNames', true, 'FileType', 'Spreadsheet', 'WriteVariableNames', false);
end