# IOC Job Analysis

Documentation for the scripts that process and visualise the results of the inverse optimal control (IOC) bilevel optimisation jobs.

## Directory layout

```
testing/job_analysis/
├── ios4/                         # IOC result .mat files — Subject 1 (Patient 4)
├── ios5/                         # IOC result .mat files — Subject 2 (Patient 5)
├── load_job_data.m               # Utility: loads all .mat files from a job directory
├── plot_main_job_graphs_s4.m     # Graphs for Subject 1 only
├── plot_main_job_graphs_s5.m     # Graphs for Subject 2 only
├── plot_main_job_graphs_s4s5.m   # Combined graphs for both subjects
├── plot_main_job_graphs_s4s5_isb_congress.m  # ISB congress variant
├── export_ioc_weights_csv.m      # Exports IOC weights to CSV (all 24 row orderings)
└── plotQuantilesAndPredictions.m # Helper for force-vs-time figures
```

## Raw IOC result files

Each `.mat` file in `ios4/` and `ios5/` stores the result for one (phase, speed, leg) condition.
Filename convention: `{phase}-{speed}-{leg}.mat`

| Token   | Values | Meaning                                  |
|---------|--------|------------------------------------------|
| `phase` | 1, 2   | 1 = stance phase (samples 1–61), 2 = swing phase (samples 61–101) |
| `speed` | 1–5    | Speed index; see speed values below      |
| `leg`   | 1, 2   | 1 = non-paretic leg, 2 = paretic leg     |

Speed values by subject:

| Speed index | Subject 1 (Patient 4) | Subject 2 (Patient 5) |
|:-----------:|:---------------------:|:---------------------:|
| 1           | 0.40 m/s              | 0.25 m/s              |
| 2           | 0.50 m/s              | 0.35 m/s              |
| 3           | 0.60 m/s              | 0.45 m/s              |
| 4           | 0.70 m/s              | 0.55 m/s              |
| 5           | 0.80 m/s              | 0.65 m/s              |

Each `.mat` file contains a struct with the following fields:

| Field   | Size     | Description                                              |
|---------|----------|----------------------------------------------------------|
| `alpha` | [15 × 1] | IOC weight vector (θ₁ … θ₁₅), one weight per cost function |
| `err`   | scalar   | RMSE between the predicted and measured force trajectories (N) |

These files are loaded by `load_job_data(dirname)`, which returns a cell array indexed as
`trials{leg, speed, phase}`. The helper `extract_properties_from_structs(trials, 'alpha')` then
assembles all weights into a single array of size `[2, 5, 2, 15, 1]` (leg × speed × phase × weight × 1).

## Graphing scripts

All four scripts must be run from the `testing/job_analysis/` directory. They load patient data
from `../../Optimization Model Data/` and IOC results from `ios4/` or `ios5/`.

### `plot_main_job_graphs_s4.m` and `plot_main_job_graphs_s5.m`

Generate per-subject figures saved to `bilevel_optim_results/job_ioc_results/patient4/` and
`patient5/` respectively:

| Output file                        | Contents                                                  |
|------------------------------------|-----------------------------------------------------------|
| `theta-vs-speed-leg-{l}-phase-{p}` | Bar chart of IOC weights across speeds, with RMSE overlay |
| `rmse-vs-speed`                    | RMSE vs. walking speed for all leg/phase combinations     |
| `force-vs-time-leg-{l}`            | Predicted vs. measured force trajectories (best-fit trial)|

### `plot_main_job_graphs_s4s5.m`

Generates a combined 2 × 4 subplot figure for both subjects at one representative speed per
subject (speed index 2 for S1, speed index 3 for S2). Saved as:

- `bilevel_optim_results/job_ioc_results/theta-vs-phase-leg-subj.pdf/png`

### `plot_main_job_graphs_s4s5_isb_congress.m`

A variant of the above formatted for ISB congress presentation. Saved as:

- `bilevel_optim_results/job_ioc_results/theta-vs-phase-leg-subj-isb-congress.pdf/png`

## Exporting IOC weights to CSV

`export_ioc_weights_csv.m` exports the full set of IOC weights and RMSEs to tabular CSV files.
Run it from the `testing/job_analysis/` directory.

Output: 24 CSV files in `bilevel_optim_results/job_ioc_results/csv_weights/`, one per ordering
permutation. Each file has 41 rows × 17 columns:

- **Row 1**: column headers  
- **Rows 2–41**: one row per condition (2 subjects × 2 legs × 2 phases × 5 speeds = 40 conditions)

| Column        | Content                                        |
|---------------|------------------------------------------------|
| 1 `Condition` | Label: `S{n}_{leg}_{phase}_{speed}m/s`         |
| 2–16 `theta_k`| IOC weight θ_k for cost function k (k = 1…15) |
| 17 `RMSE`     | Fit RMSE in Newtons                            |

The 24 files differ only in the row ordering. Each filename encodes the nesting order used to
sort the rows, e.g. `ioc_weights_subj-leg-phase-speed.csv` groups rows first by subject, then
by leg, then by phase, with speed varying fastest.
