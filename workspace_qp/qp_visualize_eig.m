function varargout = qp_visualize_eig(Q, varargin)
%QP_VISUALIZE_EIG  Heatmap of the eigendecomposition of the QP cost matrix Q.
%
%   qp_visualize_eig(Q)
%   qp_visualize_eig(Q, labels)
%   qp_visualize_eig(Q, labels, precision)
%   qp_visualize_eig(Q, labels, precision, min_cell_px)
%   qp_visualize_eig(Q, labels, precision, min_cell_px, n_levels)
%   fig = qp_visualize_eig(...)
%
%   Computes the eigendecomposition Q = V * diag(lambda) * V' and plots the
%   n×n matrix whose i-th row is  sqrt(lambda_i) * v_i, sorted by descending
%   eigenvalue.  Scaling by sqrt(lambda_i) encodes both direction and cost
%   magnitude: rows corresponding to near-zero eigenvalues are automatically
%   suppressed, while dominant modes stand out.  The colormap is symmetric
%   around zero so that positive and negative loadings are visually distinct.
%
%   The figure is sized so that every cell is at least min_cell_px × min_cell_px
%   pixels.  'PaperPositionMode' = 'auto' ensures correct export via
%   exportgraphics / saveas.
%
%   Inputs:
%     Q           - [n x n] symmetric positive-semi-definite cost matrix
%     labels      - [1 x n] string array of variable names (default: "1".."n")
%     precision   - integer, decimal places in cell text              (default: 2)
%     min_cell_px - minimum cell edge length in pixels                (default: 50)
%     n_levels    - quantization levels before plotting; [] = disabled (default: [])
%
%   Output:
%     fig         - handle to the created figure (optional)

% ---- Input validation ---------------------------------------------------
if ~ismatrix(Q) || size(Q,1) ~= size(Q,2)
    error('qp_visualize_eig:InvalidInput', 'Q must be a square matrix.');
end
n = size(Q, 1);

% ---- Optional arguments -------------------------------------------------
labels      = string(1:n);
precision   = 2;
min_cell_px = 50;
n_levels    = [];
if nargin >= 2 && ~isempty(varargin{1}), labels      = string(varargin{1}); end
if nargin >= 3 && ~isempty(varargin{2}), precision   = varargin{2};         end
if nargin >= 4 && ~isempty(varargin{3}), min_cell_px = varargin{3};         end
if nargin >= 5 && ~isempty(varargin{4}), n_levels    = varargin{4};         end

fmt = sprintf('%%.%df', precision);

% ---- Eigendecomposition, sorted descending by eigenvalue ----------------
% Symmetrize before decomposing to guarantee real output despite floating-
% point asymmetry in Q (e.g. from L*L' where L was reconstructed from a
% parameter vector).
[V, D]        = eig((Q + Q') / 2);
lambda        = real(diag(D));       % discard any residual imaginary noise
V             = real(V);
[lambda, ord] = sort(lambda, 'descend');
V             = V(:, ord);

% ---- Build scaled matrix: row i = sqrt(lambda_i) * v_i ------------------
% Clamp small negative eigenvalues (numerical noise) to zero before sqrt.
lambda_pos  = max(lambda, 0);
scaled      = diag(sqrt(lambda_pos)) * V';   % [n x n], row i = sqrt(lambda_i)*v_i

% ---- Optional quantization ----------------------------------------------
if ~isempty(n_levels) && n_levels > 1
    lo       = min(scaled(:));
    hi       = max(scaled(:));
    levels   = linspace(lo, hi, n_levels);
    [~, idx] = min(abs(scaled(:) - levels), [], 2);
    scaled   = reshape(levels(idx), n, n);
end

% ---- Figure sized so every cell is at least min_cell_px x min_cell_px --
fig  = figure('Units', 'pixels', 'PaperPositionMode', 'auto');
pos  = get(fig, 'Position');
pos(3) = n * min_cell_px;
pos(4) = n * min_cell_px;
set(fig, 'Position', pos);

ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.05 0.05 0.85 0.9]);

% ---- Heatmap with symmetric color limits --------------------------------
imagesc(ax, scaled);
abs_max  = max(abs(scaled(:)));
if abs_max > 0
    ax.CLim = [-abs_max, abs_max];
end
% Diverging blue-white-red colormap: negative loadings blue, positive red
n_cm   = 256;
r_ramp = [linspace(0, 1, n_cm/2), ones(1, n_cm/2)];
g_ramp = [linspace(0, 1, n_cm/2), linspace(1, 0, n_cm/2)];
b_ramp = [ones(1, n_cm/2),        linspace(1, 0, n_cm/2)];
colormap(ax, [r_ramp', g_ramp', b_ramp']);
colorbar(ax);
ax.YDir      = 'reverse';          % eigenvector 1 (largest lambda) at top
ax.Color     = get(fig, 'Color');
ax.TickLength = [0 0];

% ---- Axis labels --------------------------------------------------------
ax.TickLabelInterpreter = 'latex';   % set before assigning label strings
ax.XTick      = 1:n;
ax.YTick      = 1:n;
ax.XTickLabel = labels;

% Y labels: eigenvector index and its eigenvalue
y_labels      = arrayfun(@(k) sprintf('$v_{%d}$  ($\\lambda$=%.2g)', k, lambda(k)), ...
    1:n, 'UniformOutput', false);
ax.YTickLabel = y_labels;

xlabel(ax, 'Force component', 'Interpreter', 'latex');
title(ax, ['QP cost eigenmodes  ' ...
    '($\sqrt{\lambda_i}\,\mathbf{v}_i$, sorted by descending $\lambda$)'], ...
    'Interpreter', 'latex');

% ---- Cell text annotations (black on near-white, white on saturated) ----
clims      = ax.CLim;
clim_range = max(clims(2) - clims(1), eps);

for i = 1:n
    for j = 1:n
        val = scaled(i, j);
        rel = (val - clims(1)) / clim_range;   % 0=blue extreme, 1=red extreme
        % White region is near the midpoint (rel ~ 0.5); use dark text there,
        % white text near the saturated colour extremes.
        dist_from_mid = abs(rel - 0.5);
        if dist_from_mid > 0.25
            txt_color = [1 1 1];
        else
            txt_color = [0 0 0];
        end
        text(ax, j, i, sprintf(fmt, val), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'middle', ...
            'Color',               txt_color, ...
            'FontSize',            8);
    end
end

axis(ax, [0.5, n+0.5, 0.5, n+0.5]);

if nargout > 0
    varargout{1} = fig;
end

end
