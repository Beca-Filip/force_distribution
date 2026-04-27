#!/usr/bin/env python3
# export_ioc_weights_tex.py
#
# Converts each IOC-weights CSV (from export_ioc_weights_csv.m) into a
# colour-coded LaTeX table fragment.
#
# Colour rules (applied to the raw column maximum across all 40 rows):
#   grey   (gray!50)          : max < THRESH_ZERO (0.005)
#   yellow (yellow!20..!40)   : THRESH_ZERO <= max < dominant_thresholds[0]
#                               shade is determined by how many small_thresholds
#                               the column max crosses (more crossings -> darker)
#   green  (green!20..!40)    : max >= dominant_thresholds[0]
#                               shade is determined by how many dominant_thresholds
#                               the column max crosses (more crossings -> darker)
#
# Usage:
#   python3 export_ioc_weights_tex.py --small 0.15 0.30 --dominant 0.50 0.70 0.90
#
# Required LaTeX packages in the parent document:
#   \usepackage[table]{xcolor}   % columncolor, rowcolor
#   \usepackage{amsmath, bm}     % boldsymbol
#
# Run from : testing/job_analysis/
# Input    : ../../bilevel_optim_results/job_ioc_results/csv_weights/*.csv
# Output   : ../../bilevel_optim_results/job_ioc_results/tex_tables/<threshold-tag>/

import argparse
import csv
import os

INPUT_DIR_BASE  = os.path.join('..', '..', 'bilevel_optim_results',
                               'job_ioc_results', 'csv_weights')
OUTPUT_DIR_BASE = os.path.join('..', '..', 'bilevel_optim_results',
                               'job_ioc_results', 'tex_tables')

THRESH_ZERO = 0.005   # fixed: below this the column is always grey

YELLOW_LO, YELLOW_HI = 20, 40   # xcolor intensity range for yellow
GREEN_LO,  GREEN_HI  = 20, 40   # xcolor intensity range for green


# ---------------------------------------------------------------------------
# Colour logic
# ---------------------------------------------------------------------------

def _interp_intensity(level, n_levels, lo, hi):
    """
    Map an integer level in [0, n_levels] to an xcolor integer intensity in [lo, hi].

    With n_levels=0 (only one possible value) the single level gets intensity lo.
    This ensures the full range [lo, hi] is always used when n_levels >= 1.
    """
    if n_levels == 0:
        return lo
    return round(lo + (level / n_levels) * (hi - lo))


def col_color(max_val, small_thresholds, dominant_thresholds):
    """Return the xcolor spec for a weight column given its maximum value."""
    if max_val < THRESH_ZERO:
        return 'gray!50'

    if max_val >= dominant_thresholds[0]:
        # How many dominant thresholds does this column cross? (1 .. n_dom)
        n_crossed = sum(1 for t in dominant_thresholds if max_val >= t)
        n_dom     = len(dominant_thresholds)
        # level ∈ [0, n_dom-1]: 0 = crossed only the first, n_dom-1 = crossed all
        level     = n_crossed - 1
        intensity = _interp_intensity(level, n_dom - 1, GREEN_LO, GREEN_HI)
        return f'green!{intensity}'

    # Yellow band: how many small thresholds does this column cross? (0 .. n_small)
    n_crossed = sum(1 for t in small_thresholds if max_val >= t)
    n_small   = len(small_thresholds)
    intensity = _interp_intensity(n_crossed, n_small, YELLOW_LO, YELLOW_HI)
    return f'yellow!{intensity}'


# ---------------------------------------------------------------------------
# Table helpers
# ---------------------------------------------------------------------------

def make_col_spec(weight_colors):
    parts = [r'>{\columncolor{gray!10}}l']
    for c in weight_colors:
        parts.append(r'>{\columncolor{' + c + r'}}c')
    parts.append(r'>{\columncolor{gray!10}}c')          # RMSE column
    return '|' + '|'.join(parts) + '|'


def condition_to_latex(label):
    # 'S1_NP_Stance_0.40m/s' -> bold-omega superscript label in math mode
    subj, leg, phase_str, speed_str = label.split('_')
    phase = 'Stn' if phase_str == 'Stance' else 'Swg'
    speed = speed_str.replace('m/s', '')
    sup = r'\mathrm{' + subj + ',' + leg + ',' + phase + '}'
    sub = r'\mathrm{' + speed + '}'
    return r'$\boldsymbol{\omega}^{' + sup + '}_{' + sub + '}$'


def fmt_weight(val):
    r = round(float(val), 2)
    return '0' if r == 0.0 else f'{r:.2f}'


def fmt_rmse(val):
    return f'{float(val):.2f}'


def threshold_tag(small_thresholds, dominant_thresholds):
    """Short filesystem-safe tag encoding the threshold values."""
    fmt = lambda ts: '-'.join(str(int(round(t * 100))) for t in ts)
    return f's{fmt(small_thresholds)}_d{fmt(dominant_thresholds)}'


# ---------------------------------------------------------------------------
# Table builder
# ---------------------------------------------------------------------------

def build_tex(csv_path, perm_str, small_thresholds, dominant_thresholds):
    with open(csv_path, newline='') as fh:
        all_rows = list(csv.reader(fh))

    data        = all_rows[1:]
    weight_cols = [[float(row[c + 1]) for row in data] for c in range(15)]
    max_per_col = [max(col) for col in weight_cols]
    colors      = [col_color(m, small_thresholds, dominant_thresholds) for m in max_per_col]
    col_spec    = make_col_spec(colors)

    perm_display  = perm_str.replace('-', r' $\to$ ')
    small_str     = ', '.join(str(t) for t in small_thresholds)
    dominant_str  = ', '.join(str(t) for t in dominant_thresholds)

    header_cells = (
        ['~']
        + [f'$\\omega_{{{k}}}$' for k in range(1, 16)]
        + ['RMSE']
    )

    lines = []
    lines.append(r'\begin{table}[!ht]')
    lines.append(r'    \centering')
    lines.append(
        r'    \caption{IOC weight vectors $\boldsymbol{\omega}$ and fit RMSE for all'
        r' 40 conditions (2 subjects $\times$ 2 legs $\times$ 2 phases $\times$ 5 speeds).'
        r' Row ordering: ' + perm_display + r'.'
        r' Column shading: grey $=$ zero everywhere ($<' + str(THRESH_ZERO) + r'$);'
        r' yellow shades $=$ small contribution (thresholds: ' + small_str + r');'
        r' green shades $=$ dominant contribution (thresholds: ' + dominant_str + r').}'
    )
    lines.append(r'    \label{tab:ioc-weights-' + perm_str + r'}')
    lines.append(r'    \footnotesize')
    lines.append(r'    \setlength{\tabcolsep}{3.2pt}')
    lines.append(r'    \renewcommand{\arraystretch}{1.05}')
    lines.append(r'    \begin{tabular}{' + col_spec + r'}')
    lines.append(r'    \hline')
    lines.append(r'    \rowcolor{gray!10}')
    lines.append('        ' + ' & '.join(header_cells) + r' \\ \hline')

    for row in data:
        label = condition_to_latex(row[0])
        cells = (
            [label]
            + [fmt_weight(row[c + 1]) for c in range(15)]
            + [fmt_rmse(row[16])]
        )
        lines.append(r'        \rule{0pt}{13pt}')
        lines.append('        ' + ' & '.join(cells) + r' \\ \hline')

    lines.append(r'    \end{tabular}')
    lines.append(r'\end{table}')

    return '\n'.join(lines) + '\n'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert IOC-weights CSVs to LaTeX tables with colour-coded columns.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Example:\n'
            '  python3 export_ioc_weights_tex.py --small 0.15 0.30 --dominant 0.50 0.70 0.90\n\n'
            'Colour bands:\n'
            '  grey   : column max < 0.005\n'
            '  yellow : column max >= 0.005 AND < first dominant threshold\n'
            '           shade = yellow!20 (none crossed) .. yellow!40 (all small thresholds crossed)\n'
            '  green  : column max >= first dominant threshold\n'
            '           shade = green!20 (only first crossed) .. green!40 (all dominant crossed)\n'
        )
    )
    parser.add_argument(
        '--small', '-s', nargs='+', type=float, required=True, metavar='T',
        help='Upper-bound thresholds for yellow sub-bands (increasing order)'
    )
    parser.add_argument(
        '--dominant', '-d', nargs='+', type=float, required=True, metavar='T',
        help='Lower-bound thresholds for green sub-bands (increasing order)'
    )
    parser.add_argument(
        '--input-dir', '-i', default=INPUT_DIR_BASE,
        help='Directory containing the CSV files (default: %(default)s)'
    )
    parser.add_argument(
        '--output-dir', '-o', default=None,
        help='Output directory for .tex files (default: auto-named under tex_tables/)'
    )
    args = parser.parse_args()

    small_thresholds    = sorted(args.small)
    dominant_thresholds = sorted(args.dominant)

    if small_thresholds[-1] >= dominant_thresholds[0]:
        parser.error(
            f'All small thresholds must be strictly below all dominant thresholds '
            f'(largest small = {small_thresholds[-1]}, smallest dominant = {dominant_thresholds[0]}).'
        )

    output_dir = (args.output_dir or
                  os.path.join(OUTPUT_DIR_BASE, threshold_tag(small_thresholds,
                                                               dominant_thresholds)))
    os.makedirs(output_dir, exist_ok=True)

    print(f'Small thresholds    : {small_thresholds}')
    print(f'Dominant thresholds : {dominant_thresholds}')
    print(f'Output              : {output_dir}')

    # Print the expected colour mapping for transparency
    print()
    print('Yellow sub-bands:')
    bounds = [THRESH_ZERO] + small_thresholds + [dominant_thresholds[0]]
    n_small = len(small_thresholds)
    for k in range(n_small + 1):
        lo_b, hi_b = bounds[k], bounds[k + 1]
        intensity  = _interp_intensity(k, n_small, YELLOW_LO, YELLOW_HI)
        hi_str     = f'{hi_b}' if k < n_small else f'{hi_b} (exclusive)'
        print(f'  [{lo_b}, {hi_str})  ->  yellow!{intensity}')

    print('Green sub-bands:')
    n_dom   = len(dominant_thresholds)
    d_bounds = dominant_thresholds + [float('inf')]
    for k in range(n_dom):
        lo_b, hi_b = d_bounds[k], d_bounds[k + 1]
        intensity  = _interp_intensity(k, n_dom - 1, GREEN_LO, GREEN_HI)
        hi_str     = f'{hi_b}' if hi_b != float('inf') else 'inf'
        print(f'  [{lo_b}, {hi_str})  ->  green!{intensity}')
    print()

    csv_files = sorted(f for f in os.listdir(args.input_dir) if f.endswith('.csv'))
    for fname in csv_files:
        perm_str = fname[len('ioc_weights_'):-len('.csv')]
        csv_path = os.path.join(args.input_dir, fname)
        tex_path = os.path.join(output_dir, fname[:-len('.csv')] + '.tex')

        content = build_tex(csv_path, perm_str, small_thresholds, dominant_thresholds)
        with open(tex_path, 'w') as fh:
            fh.write(content)
        print(f'Written: {tex_path}')

    print(f'\nDone. {len(csv_files)} .tex files written.')


if __name__ == '__main__':
    main()
