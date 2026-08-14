function [Q_to, l_to, remap_info] = remap_qp_to_subject(Q_from, l_from, names_from, names_to, varargin)
%REMAP_QP_TO_SUBJECT  Move QP cost weights (Q, l) between two muscle sets.
%
%   [Q_to, l_to] = REMAP_QP_TO_SUBJECT(Q_from, l_from, names_from, names_to)
%   [Q_to, l_to, remap_info] = REMAP_QP_TO_SUBJECT(...)
%   [...] = REMAP_QP_TO_SUBJECT(..., 'q_fill', 'median', 'l_fill', 'zero')
%
%   Weights identified on one subject are expressed on another subject's
%   muscle set.  Matching is by NAME (case-insensitive, whitespace-trimmed),
%   so it handles both re-ordering and differing muscle counts.
%
%   For patients 4 and 5 the sets are:
%     35 muscles (P4), 34 muscles (P5), 33 in common,
%     edl + fdl exist only in P4,  tfl exists only in P5.
%
%   THREE CASES
%   -----------
%     matched  - present in both.  The entire matched sub-block of Q and the
%                matched entries of l are carried over EXACTLY, re-indexed
%                into the target ordering.
%     dropped  - present in names_from only.  Their rows/columns are removed.
%                Safe: any principal submatrix of a PSD matrix is PSD.
%     inserted - present in names_to only.  There is no identified weight for
%                them, so they receive a positive DIAGONAL entry and ZERO
%                coupling to every other muscle (see WHY DIAGONAL below).
%
%   WHY DIAGONAL-ONLY FILL
%   ----------------------
%   The obvious alternative -- copy an anatomically similar muscle's whole
%   row/column -- is unsafe.  Two identical rows make Q singular, which is
%   exactly the rank loss that motivated abandoning the synergy
%   parameterisation: a singular Q makes the QP's minimiser non-unique, so
%   the inverse problem stops being well posed.  A diagonal fill cannot do
%   that.  In the target ordering the result is similarity-equivalent to
%
%       [ Q_matched   0       ]
%       [ 0           diag(d) ]
%
%   whose eigenvalues are eig(Q_matched) together with the fill values
%   themselves.  So PSD is preserved whenever Q_from was PSD and every fill
%   value is > 0, and no new null direction is ever created.
%
%   Interpretation: an inserted muscle is penalised on its own normalised
%   deviation and is a priori uncoupled from the rest -- a neutral stance,
%   not a claim that the muscle is unimportant.
%
%   OPTIONS (name/value)
%   --------------------
%     'q_fill'  Diagonal value for inserted muscles. Default 'median'.
%               'median' | 'mean'  - statistic of diag(Q_from) over the
%                                    MATCHED muscles only (so the fill is on
%                                    the same scale as the weights kept).
%               numeric scalar     - use this value directly (must be > 0).
%               'donor'            - use diag(Q_from) of a named donor muscle
%                                    given in 'donor_map'.
%     'l_fill'  Linear-term value for inserted muscles. Default 'zero'.
%               'zero'             - no linear bias.  With zero coupling this
%                                    drives the muscle toward f_norm = 0,
%                                    i.e. toward the subject's own mean force
%                                    -- the natural "no information" default.
%               'median' | 'mean'  - statistic of l_from over matched muscles.
%               numeric scalar | 'donor'  - as above.
%     'donor_map'  {new_name, donor_name; ...} cell array, required when
%               either fill is 'donor'.  Donor must exist in names_from.
%     'psd_tol' Relative tolerance for the PSD check. Default 1e-10.
%     'verbose' Print a summary. Default false.
%
%   IMPORTANT CAVEAT -- ORDER-DEPENDENT WHITENING
%   ---------------------------------------------
%   Q does not act on raw forces.  The cost is
%       0.5*f_norm'*Q*f_norm + l'*f_norm,   f_norm = W*(f - mu),
%   and compute_F_invcov builds W = inv(chol(Sigma,'lower')), which is LOWER
%   TRIANGULAR and therefore ORDER-DEPENDENT: coordinate i of f_norm is
%   muscle i's residual after regressing out muscles 1..i-1.  Cholesky is not
%   permutation-equivariant, so "coordinate i" does not mean the same thing
%   for two subjects with different muscle orderings.  Inserting tfl at P5
%   index 18 shifts the conditioning set of every later coordinate.
%
%   This function remaps INDICES correctly and exactly.  It cannot repair a
%   mismatch in what those indices mean.  Treat a cross-subject transfer
%   under Cholesky whitening as an approximation.  A symmetric whitening
%   (Sigma^-1/2) or a diagonal one (per-muscle z-score) IS permutation-
%   equivariant and would make this remap exact; both require changing
%   compute_F_invcov and re-running the fits.
%
%   SCALE
%   -----
%   The QP minimiser is invariant to (Q, l) -> (c*Q, c*l) for any c > 0, so
%   the weights are identified only up to positive scale.  Every statistic
%   fill option here is homogeneous of degree 1, so this function commutes
%   with rescaling: remap(c*Q, c*l) == c*remap(Q, l).
%
%   See also CHOL_VEC_TO_Q, COMPUTE_F_INVCOV, QP_SUBROUTINE.

% ---- Parse options ------------------------------------------------------
p = inputParser;
p.FunctionName = 'remap_qp_to_subject';
addParameter(p, 'q_fill',    'median');
addParameter(p, 'l_fill',    'zero');
addParameter(p, 'donor_map', {});
addParameter(p, 'psd_tol',   1e-10, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'verbose',   false);
parse(p, varargin{:});

q_fill    = p.Results.q_fill;
l_fill    = p.Results.l_fill;
donor_map = p.Results.donor_map;
psd_tol   = p.Results.psd_tol;
verbose   = logical(p.Results.verbose);

% ---- Normalise and validate the name lists ------------------------------
names_from = local_to_cellstr(names_from, 'names_from');
names_to   = local_to_cellstr(names_to,   'names_to');

n_from = numel(names_from);
n_to   = numel(names_to);

if ~isequal(size(Q_from), [n_from, n_from])
    error('remap_qp_to_subject:sizeMismatch', ...
        'Q_from is %dx%d but names_from has %d entries.', ...
        size(Q_from, 1), size(Q_from, 2), n_from);
end
l_from = l_from(:);
if numel(l_from) ~= n_from
    error('remap_qp_to_subject:sizeMismatch', ...
        'l_from has %d entries but names_from has %d.', numel(l_from), n_from);
end

key_from = lower(strtrim(names_from));
key_to   = lower(strtrim(names_to));

local_assert_unique(key_from, 'names_from');
local_assert_unique(key_to,   'names_to');

% ---- Match by name ------------------------------------------------------
[is_matched, loc_in_from] = ismember(key_to, key_from);

idx_to_matched    = find(is_matched);
idx_from_matched  = loc_in_from(is_matched);
idx_to_inserted   = find(~is_matched);
idx_from_dropped  = setdiff((1:n_from)', idx_from_matched(:));

n_matched  = numel(idx_to_matched);
n_inserted = numel(idx_to_inserted);
n_dropped  = numel(idx_from_dropped);

if n_matched == 0
    error('remap_qp_to_subject:noOverlap', ...
        'No muscle names are common to the two sets -- nothing to transfer.');
end

% ---- Fill values --------------------------------------------------------
% Statistics are taken over the MATCHED muscles only, so the fill sits on the
% same scale as the weights actually being carried across.
diag_matched = diag(Q_from(idx_from_matched, idx_from_matched));
l_matched    = l_from(idx_from_matched);

q_fill_values = local_resolve_fill(q_fill, 'q_fill', n_inserted, ...
    diag_matched, names_to(idx_to_inserted), key_from, ...
    diag(Q_from), donor_map);

l_fill_values = local_resolve_fill(l_fill, 'l_fill', n_inserted, ...
    l_matched, names_to(idx_to_inserted), key_from, ...
    l_from, donor_map);

if any(q_fill_values <= 0)
    error('remap_qp_to_subject:nonPositiveFill', ...
        ['q_fill produced a non-positive diagonal value (%g). A positive ' ...
         'value is required to keep Q positive definite on the inserted ' ...
         'muscles.'], min(q_fill_values));
end

% ---- Build the remapped weights ----------------------------------------
Q_to = zeros(n_to, n_to);
Q_to(idx_to_matched, idx_to_matched) = Q_from(idx_from_matched, idx_from_matched);

% Inserted muscles: positive diagonal, zero coupling (rows/cols already 0).
for ii = 1:n_inserted
    Q_to(idx_to_inserted(ii), idx_to_inserted(ii)) = q_fill_values(ii);
end

% Kill any asymmetry inherited from a slightly non-symmetric Q_from.
Q_to = (Q_to + Q_to') / 2;

l_to = zeros(n_to, 1);
l_to(idx_to_matched)  = l_from(idx_from_matched);
l_to(idx_to_inserted) = l_fill_values;

% ---- Diagnostics --------------------------------------------------------
eig_from = eig((Q_from + Q_from') / 2);
eig_to   = eig(Q_to);

remap_info = struct();
remap_info.names_from       = names_from;
remap_info.names_to         = names_to;
remap_info.names_matched    = names_to(idx_to_matched);
remap_info.names_dropped    = names_from(idx_from_dropped);
remap_info.names_inserted   = names_to(idx_to_inserted);
remap_info.idx_from_matched = idx_from_matched(:);
remap_info.idx_to_matched   = idx_to_matched(:);
remap_info.idx_from_dropped = idx_from_dropped(:);
remap_info.idx_to_inserted  = idx_to_inserted(:);
remap_info.n_matched        = n_matched;
remap_info.n_dropped        = n_dropped;
remap_info.n_inserted       = n_inserted;
remap_info.q_fill_values    = q_fill_values(:);
remap_info.l_fill_values    = l_fill_values(:);
remap_info.eig_min_from     = min(eig_from);
remap_info.eig_min_to       = min(eig_to);
remap_info.cond_from        = max(eig_from) / max(min(eig_from), eps);
remap_info.cond_to          = max(eig_to)   / max(min(eig_to),   eps);

% PSD check on the RESULT.  A negative eigenvalue here means Q_from was not
% PSD to begin with (the construction cannot introduce one).
scale_to = max(abs(eig_to));
remap_info.is_psd = min(eig_to) >= -psd_tol * max(scale_to, 1);

if ~remap_info.is_psd
    warning('remap_qp_to_subject:notPSD', ...
        ['Remapped Q is not PSD (min eigenvalue %.3e). The matched block of ' ...
         'Q_from must already have been indefinite -- this construction ' ...
         'cannot create a negative eigenvalue.'], min(eig_to));
end

if verbose
    fprintf('remap_qp_to_subject: %d -> %d muscles\n', n_from, n_to);
    fprintf('  matched  : %d\n', n_matched);
    fprintf('  dropped  : %d%s\n', n_dropped, local_name_list(remap_info.names_dropped));
    fprintf('  inserted : %d%s\n', n_inserted, local_name_list(remap_info.names_inserted));
    for ii = 1:n_inserted
        fprintf('      %-12s Q(i,i) = %.4g,  l(i) = %.4g\n', ...
            remap_info.names_inserted{ii}, q_fill_values(ii), l_fill_values(ii));
    end
    fprintf('  min eig  : %.3e -> %.3e\n', remap_info.eig_min_from, remap_info.eig_min_to);
    fprintf('  cond(Q)  : %.3e -> %.3e\n', remap_info.cond_from, remap_info.cond_to);
end

end

% =========================================================================
% Helpers
% =========================================================================

function names = local_to_cellstr(names, arg_name)
%LOCAL_TO_CELLSTR  Accept cellstr, string array or char matrix.
if isstring(names) || ischar(names)
    names = cellstr(names);
end
if ~iscellstr(names) %#ok<ISCLSTR>
    error('remap_qp_to_subject:badNames', ...
        '%s must be a cell array of char vectors, a string array or a char matrix.', ...
        arg_name);
end
names = names(:);
end

% -------------------------------------------------------------------------

function local_assert_unique(keys, arg_name)
%LOCAL_ASSERT_UNIQUE  Duplicate names would make the match ambiguous.
[unique_keys, ~, group] = unique(keys);
counts = accumarray(group, 1);
if any(counts > 1)
    dupes = unique_keys(counts > 1);
    error('remap_qp_to_subject:duplicateNames', ...
        '%s contains duplicate muscle name(s): %s', ...
        arg_name, strjoin(dupes(:)', ', '));
end
end

% -------------------------------------------------------------------------

function values = local_resolve_fill(spec, spec_name, n_inserted, ...
    matched_values, inserted_names, key_from, from_values, donor_map)
%LOCAL_RESOLVE_FILL  Turn a fill specification into one value per inserted muscle.

if n_inserted == 0
    values = zeros(0, 1);
    return
end

if isnumeric(spec)
    if ~isscalar(spec)
        error('remap_qp_to_subject:badFill', ...
            'Numeric %s must be a scalar.', spec_name);
    end
    values = repmat(spec, n_inserted, 1);
    return
end

if ~(ischar(spec) || isstring(spec))
    error('remap_qp_to_subject:badFill', ...
        '%s must be a numeric scalar or one of ''median'', ''mean'', ''zero'', ''donor''.', ...
        spec_name);
end

switch lower(char(spec))
    case 'median'
        values = repmat(median(matched_values), n_inserted, 1);
    case 'mean'
        values = repmat(mean(matched_values), n_inserted, 1);
    case 'zero'
        values = zeros(n_inserted, 1);
    case 'donor'
        values = zeros(n_inserted, 1);
        if isempty(donor_map)
            error('remap_qp_to_subject:missingDonorMap', ...
                '%s = ''donor'' requires a non-empty ''donor_map''.', spec_name);
        end
        donor_keys  = lower(strtrim(donor_map(:, 1)));
        donor_names = lower(strtrim(donor_map(:, 2)));
        for ii = 1:n_inserted
            this_key = lower(strtrim(inserted_names{ii}));
            row = find(strcmp(donor_keys, this_key), 1);
            if isempty(row)
                error('remap_qp_to_subject:missingDonor', ...
                    'donor_map has no entry for inserted muscle ''%s''.', ...
                    inserted_names{ii});
            end
            donor_idx = find(strcmp(key_from, donor_names{row}), 1);
            if isempty(donor_idx)
                error('remap_qp_to_subject:unknownDonor', ...
                    'Donor muscle ''%s'' (for ''%s'') is not in names_from.', ...
                    donor_map{row, 2}, inserted_names{ii});
            end
            values(ii) = from_values(donor_idx);
        end
    otherwise
        error('remap_qp_to_subject:badFill', ...
            'Unknown %s option ''%s''.', spec_name, char(spec));
end
end

% -------------------------------------------------------------------------

function s = local_name_list(names)
if isempty(names)
    s = '';
else
    s = sprintf(' (%s)', strjoin(names(:)', ', '));
end
end
