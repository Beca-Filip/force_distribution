function varargout = qp_visualize_weights(Q, l, varargin)
%QP_VISUALIZE_WEIGHTS  Heatmap of QP cost weights Q (quadratic) and l (linear).
%
%   qp_visualize_weights(Q, l)
%   qp_visualize_weights(Q, l, labels)
%   qp_visualize_weights(Q, l, labels, precision)
%   qp_visualize_weights(Q, l, labels, precision, minCellPx)
%   fig = qp_visualize_weights(...)
%
%   Renders a single merged heatmap with n+1 columns and n rows:
%     Column 1       : linear weight vector l (one value per row/feature)
%     Columns 2:n+1  : lower triangle of Q; upper triangle is left blank (NaN)
%
%   A vertical white separator line accentuates the boundary between l and Q.
%   Each visible cell carries a text annotation of its numeric value.
%   The colorbar is shared across both blocks.
%
%   The figure is sized so that every cell is at least minCellPx × minCellPx
%   pixels, allowing off-screen rendering larger than the physical monitor.
%   'PaperPositionMode' = 'auto' ensures exportgraphics / saveas capture the
%   full canvas without clipping.
%
%   Inputs:
%     Q           - [n x n] symmetric positive-semi-definite weight matrix
%     l           - [n x 1] or [1 x n] linear weight vector
%     labels      - [1 x n] string array of variable names (default: "1".."n")
%     precision   - integer, decimal places shown in each cell (default: 2)
%     min_cell_px - minimum cell edge length in pixels            (default: 50)
%     n_levels    - number of quantization levels applied to all values before
%                   plotting; [] or 0 disables quantization               (default: [])
%                   Example: n_levels=3 maps every value to the nearest of
%                   {min_val, mid_val, max_val} across all of l and Q.
%
%   Output:
%     fig         - handle to the created figure (optional)

% ---- Input validation ---------------------------------------------------
if ~ismatrix(Q) || size(Q,1) ~= size(Q,2)
    error('qp_visualize_weights:InvalidInput', 'Q must be a square matrix.');
end
n = size(Q, 1);
if ~isvector(l) || numel(l) ~= n
    error('qp_visualize_weights:InvalidInput', ...
        'l must be a vector with the same number of elements as Q.');
end

% ---- Optional arguments -------------------------------------------------
labels      = string(1:n);
precision   = 2;
min_cell_px = 50;
n_levels    = [];
if nargin >= 3 && ~isempty(varargin{1}), labels      = string(varargin{1}); end
if nargin >= 4 && ~isempty(varargin{2}), precision   = varargin{2};         end
if nargin >= 5 && ~isempty(varargin{3}), min_cell_px = varargin{3};         end
if nargin >= 6 && ~isempty(varargin{4}), n_levels    = varargin{4};         end

fmt = sprintf('%%.%df', precision);

% ---- Build combined matrix (n rows x n+1 cols) --------------------------
% Transpose of [l; Q] gives:
%   col 1   = l  (one entry per feature row)
%   col j+1 = j-th column of Q  (Q is symmetric so Q^T = Q)
l        = reshape(l, 1, []);   % ensure row vector before stacking
combined = [l; Q].';            % (n+1 x n).' -> n x (n+1)

% Blank the strict upper triangle of Q inside the combined matrix.
% combined(i, j+1) = Q(j,i); blanking upper-triangle of Q (j < i) means
% blanking combined(i, j+1) where i > j, i.e. lower-triangle of the Q block.
q_block = combined(:, 2:end);
q_block(tril(true(n), -1)) = NaN;
combined(:, 2:end) = q_block;

% ---- Optional quantization of visible values ----------------------------
if ~isempty(n_levels) && n_levels > 1
    valid_mask = ~isnan(combined);
    vals       = combined(valid_mask);          % column vector of valid values
    lo         = min(vals);
    hi         = max(vals);
    levels     = linspace(lo, hi, n_levels);    % [1 x n_levels] quantization grid
    % Map each value to the nearest level via distance matrix [numel(vals) x n_levels]
    [~, idx]   = min(abs(vals - levels), [], 2);
    combined(valid_mask) = levels(idx);
end

% ---- Figure sized so every cell is at least min_cell_px x min_cell_px --
n_rows = n;
n_cols = n + 1;
fig    = figure('Units', 'pixels', 'PaperPositionMode', 'auto');
pos    = get(fig, 'Position');
pos(3) = n_cols * min_cell_px;
pos(4) = n_rows * min_cell_px;
set(fig, 'Position', pos);

ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);

% ---- Heatmap ------------------------------------------------------------
imagesc(ax, combined);
colormap(ax, parula);
colorbar(ax);
ax.YDir    = 'reverse';          % row 1 at the top
ax.Color   = get(fig, 'Color'); % NaN cells inherit background (white)
ax.TickLength = [0 0];

% ---- Separator between l column and Q block -----------------------------
hold(ax, 'on');
plot(ax, [1.5, 1.5], [0.5, n_rows+0.5], 'w-', 'LineWidth', 2);
hold(ax, 'off');

% ---- Axis labels --------------------------------------------------------
ax.XTick      = 1:n_cols;
ax.YTick      = 1:n_rows;
ax.XTickLabel = ["l", labels(:)'];
ax.YTickLabel = labels;

xlabel(ax, 'Terms  ({\boldmath$l$} at left, {\boldmath$Q$} columns to right)', 'Interpreter', 'latex');
title(ax,  'QP cost weights: linear {\boldmath$l$} (left) and quadratic {\boldmath$Q$} lower triangle (right)', ...
    'Interpreter', 'latex');

% ---- Cell text annotations (black on bright, white on dark) -------------
clims      = ax.CLim;
clim_range = max(clims(2) - clims(1), eps);

for i = 1:n_rows
    for j = 1:n_cols
        val = combined(i, j);
        if isnan(val), continue; end
        rel = (val - clims(1)) / clim_range;
        if rel > 0.5
            txtColor = [1 1 1];   % white on dark background
        else
            txtColor = [0 0 0];   % black on bright background
        end
        text(ax, j, i, sprintf(fmt, val), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'middle', ...
            'Color',               txtColor, ...
            'FontSize',            10);
    end
end

axis(ax, [0.5, n_cols+0.5, 0.5, n_rows+0.5]);

if nargout > 0
    varargout{1} = fig;
end

end
