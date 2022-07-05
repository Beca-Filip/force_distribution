function [h_ax, h] = patch_vector_quantities_opts_shape(t, f, h_ax, opts, pcol, shape, varargin)
%PATCH_VECTOR_QUANTITIES_OPTS automatically determines the size of the subgrid
%in which to plot a vector trajectory, and performs the plotting, while
%allowing some freedom with options.
%
%   
%   Takes in an options structure for axes parameters (title, xlabel,
%   ylabel, zlabel) where fields are named with the function they're
%   supposed to call and the field value is a function which takes in the
%   order of the trajectory to plot and returns a cell array with function
%   arguments
%   e.g. opts.title = @(n) {sprintf('trajectory %d', n)}
%
%   Takes in plot options which are common to all plots.
%   e.g.
%   [h_ax, h] = PLOT_VECTOR_QUANTITIES(t, f, h_ax, opts, 'LineWidth', 2);

% Size of quantities to plot 
% n - number of quantities
% N - number of samples per quantity
[n, N] = size(f);

% If shape is given
if ~isempty(shape)
    
    % Check if its incorrect
    if ~isnumeric(shape) || numel(shape)~=2 || any(shape <= 0) || any(shape ~= round(shape))
        error('shape parameter must be a two-element vector containing positive integers.');
    end
    % Check if its not big enough
    if prod(shape) < n
        error('shape must be big enough to contain all vector trajectories.');
    end
    
    
    % Set the number of subgrid rows and columns
    figrows = shape(1);
    figcols = shape(2);
    
% Otherwise automatically determine subgrid size
else

    % Round
    figrows = ceil(sqrt(n));
    figcols = ceil(sqrt(n));

    % Try to reduce number of columns
    if (figcols-1) * figrows >= n
        figcols = figcols - 1;
    end
end

% Create an array of graphical objects
h = gobjects(figrows, figcols);
if ~isempty(h_ax)
    h_ax = cell(figrows, figcols);
end

% Treat optional arguments
if ~isempty(opts)
    optnames = fieldnames(opts);
else
    optnames = [];
end
% Number of optional arguments
noptnames = length(optnames);

% For each row and column
for ii = 1 : figrows
    for jj = 1 : figcols
        
        % Determine the current order
        curr = (ii-1)*figcols + jj;
        % If current order surpasses the number of elements to plot
        if curr > n
            break
        end
        
        % Create current subplot
        h_ax{ii, jj} = subplot(figrows, figcols, curr);
        % Hold
        hold on;
        
        % If arguments are passed plot with them
        if ~isempty(varargin)
            
            h(ii, jj) = patch(t, f(curr, :), pcol, varargin{:});
            
        % If not dont
        else
           
            h(ii, jj) = patch(t, f(curr, :), pcol);
            
        end
        
        % If optional arguments are passed
        if ~isempty(opts)
           
            % For all passed field names
            for nopt = 1 : noptnames
                % Get arguments of the optional function in cell array form
                % as described in description
                optargs = opts.(optnames{nopt})(curr);
                
                % Check its cell array character
                if ~iscell(optargs)
                    errmsg = sprintf('opts.%s must return a cell array.', optnames{nopt});
                    error(errmsg);
                end
                
                % Evaluate the function with the function name being the
                % opts structure fieldname and the calculated arguments.
                eval(sprintf('%s(optargs{:})', optnames{nopt}));
            end
            
        end
    end
end


end