function [h, varargout] = plot_num_total_gridpoints(n_simplex, m_range, varargin)
% PLOT_NUM_TOTAL_GRIDPOINTS plots the total number of grid points of a
% simplex with a fixed amount of vertices, given a variable sized partition
% along the edges.
%   
%   h = PLOT_NUM_TOTAL_GRIDPOINTS(n_simplex, m_range)
%   [h, v] = PLOT_NUM_TOTAL_GRIDPOINTS(n_simplex, m_range)
%   h = PLOT_NUM_TOTAL_GRIDPOINTS(n_simplex, m_range, Name, Value)
%   [h, v] = PLOT_NUM_TOTAL_GRIDPOINTS(n_simplex, m_range, Name, Value)
%   Inputs:
%   n_simplex ~ number of vertices of the simplex
%   m_range ~ Range of simplex partition rank
%   Name, Value ~ arguments for the plot function
%   Outputs:
%   h ~ plot handle
%   v ~ vector of plotted values

% Round partition 
m_range = round(m_range);

% Vector of partition ranks
pr = m_range(1) : m_range(2);

% Vector of total number of gridpoints
v = zeros(range(m_range)+1, 1);

% For all ranks calculate number of grid points
for ii = 1 : length(pr)
    v(ii) = nchoosek(n_simplex + pr(ii) - 2, n_simplex - 1);
end

% If arguments are passed
if ~isempty(varargin)
    h = plot(pr, log10(v), 'Marker', 'o', varargin{:});
else
    h = plot(pr, log10(v), 'Marker', 'o');
end

% Show exact number
for ii = 1 : length(pr)
   text(pr(ii), log10(v(ii)), sprintf('%d', v(ii))); 
end

% Some aesthetics
xlabel("Partition Rank");
ylabel("# Grid Points [log_{10}]");
title({sprintf("Number of grid points of a"); sprintf(" %d-dimensional simplex for different partition ranks", n_simplex)});
grid;

% If arguments are requested
if nargout > 1
    varargout{2} = v;
end

end