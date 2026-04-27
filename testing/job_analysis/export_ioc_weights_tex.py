#!/usr/bin/env python3
# export_ioc_weights_tex.py
#
# Converts each IOC-weights CSV produced by export_ioc_weights_csv.m into a
# self-contained LaTeX table fragment with colour-coded weight columns.
#
# Colour rules (evaluated on raw, un-rounded column maxima across all 40 rows):
#   grey   (gray!50)   : max < 0.005  -- rounds to 0.00 in every condition
#   yellow (yellow!20) : 0.005 <= max <= 0.200  -- small but non-zero contribution
#   green  (green!10)  : max > 0.200  -- dominant in at least one condition
#
# Required LaTeX packages in the parent document:
#   \usepackage[table]{xcolor}   % columncolor, rowcolor
#   \usepackage{amsmath, bm}     % boldsymbol
#
# Run from : testing/job_analysis/
# Input    : ../../bilevel_optim_results/job_ioc_results/csv_weights/*.csv
# Output   : ../../bilevel_optim_results/job_ioc_results/tex_tables/*.tex

import csv
import os

INPUT_DIR  = os.path.join('..', '..', 'bilevel_optim_results',
                          'job_ioc_results', 'csv_weights')
OUTPUT_DIR = os.path.join('..', '..', 'bilevel_optim_results',
                          'job_ioc_results', 'tex_tables')

THRESH_ZERO   = 0.005   # below this -> grey
THRESH_YELLOW = 0.200   # at or below this -> yellow; above -> green


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def col_color(max_val):
    if max_val < THRESH_ZERO:
        return 'gray!50'
    if max_val <= THRESH_YELLOW:
        return 'yellow!20'
    return 'green!10'


def make_col_spec(weight_colors):
    """Build the full |col|col|...| tabular specification string."""
    parts = [r'>{\columncolor{gray!10}}l']          # condition-label column
    for c in weight_colors:
        parts.append(r'>{\columncolor{' + c + r'}}c')
    parts.append(r'>{\columncolor{gray!10}}c')      # RMSE column
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
    """Round to 2 d.p.; display bare '0' for exact zero."""
    r = round(float(val), 2)
    return '0' if r == 0.0 else f'{r:.2f}'


def fmt_rmse(val):
    return f'{float(val):.2f}'


# ---------------------------------------------------------------------------
# Table builder
# ---------------------------------------------------------------------------

def build_tex(csv_path, perm_str):
    with open(csv_path, newline='') as fh:
        all_rows = list(csv.reader(fh))

    data = all_rows[1:]   # 40 condition rows, skip header

    # Per-column statistics (theta_1 ... theta_15 are columns 1-15)
    weight_cols = [[float(row[c + 1]) for row in data] for c in range(15)]
    max_per_col = [max(col) for col in weight_cols]
    colors      = [col_color(m) for m in max_per_col]
    col_spec    = make_col_spec(colors)

    perm_display = perm_str.replace('-', r' $\to$ ')

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
        r' 40 conditions (2 subjects $\times$ 2 legs $\times$ 2 phases'
        r' $\times$ 5 speeds). Row ordering: ' + perm_display + r'.'
        r' Column shading: grey $=$ zero everywhere ($<0.005$),'
        r' yellow $=$ small peak contribution (${\leq}0.20$),'
        r' green $=$ dominant contribution ($>0.20$ for at least one condition).}'
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
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    csv_files = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith('.csv'))
    for fname in csv_files:
        perm_str = fname[len('ioc_weights_'):-len('.csv')]
        csv_path = os.path.join(INPUT_DIR, fname)
        tex_path = os.path.join(OUTPUT_DIR, fname[:-len('.csv')] + '.tex')

        content = build_tex(csv_path, perm_str)
        with open(tex_path, 'w') as fh:
            fh.write(content)
        print(f'Written: {tex_path}')

    print(f'\nDone. {len(csv_files)} .tex files written to {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
