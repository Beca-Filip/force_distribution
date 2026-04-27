% export_ioc_weights_csv.m
% Exports IOC weights and RMSEs for both subjects to CSV files.
% Generates all 24 orderings obtained by permuting the nesting of
% {subject, leg, phase, speed} when enumerating the 40 conditions.
%
% Run from: testing/job_analysis/
% Output:   ../../bilevel_optim_results/job_ioc_results/csv_weights/

close all; clear all; clc;

addpath(fullfile('..', '..', 'utils'));

%% Load IOC data
trials_s1 = load_job_data('ios4');
alpha_s1   = extract_properties_from_structs(trials_s1, 'alpha');
rmse_s1    = extract_properties_from_structs(trials_s1, 'err');

trials_s2 = load_job_data('ios5');
alpha_s2   = extract_properties_from_structs(trials_s2, 'alpha');
rmse_s2    = extract_properties_from_structs(trials_s2, 'err');

%% Speed values (m/s) for each subject
speeds_s1 = linspace(0.40, 0.80, 5);
speeds_s2 = linspace(0.25, 0.65, 5);

%% Build flat row array — canonical order: subj > leg > phase > speed
% flat index: (subj-1)*20 + (leg-1)*10 + (phase-1)*5 + speed
subj_labels  = {'S1',     'S2'  };
leg_labels   = {'NP',     'P'   };
phase_labels = {'Stance', 'Swing'};

n_rows = 40;
rows(n_rows) = struct('label', '', 'alpha', [], 'rmse', []);

row_idx = 0;
for subj = 1:2
    if subj == 1
        alpha_data = alpha_s1;
        rmse_data  = rmse_s1;
        speeds     = speeds_s1;
    else
        alpha_data = alpha_s2;
        rmse_data  = rmse_s2;
        speeds     = speeds_s2;
    end
    for leg = 1:2
        for phase = 1:2
            for speed = 1:5
                row_idx = row_idx + 1;
                rows(row_idx).label = sprintf('%s_%s_%s_%.2fm/s', ...
                    subj_labels{subj}, leg_labels{leg}, phase_labels{phase}, speeds(speed));
                rows(row_idx).alpha = reshape(alpha_data(leg, speed, phase, :, :), [], 1);
                rows(row_idx).rmse  = rmse_data(leg, speed, phase);
            end
        end
    end
end

%% Column headers
weight_col_names = arrayfun(@(k) sprintf('theta_%d', k), 1:15, 'UniformOutput', false);
col_headers = [{'Condition'}, weight_col_names, {'RMSE'}];

%% Generate and write all 24 permutations
factor_names  = {'subj', 'leg', 'phase', 'speed'};
factor_ranges = {1:2, 1:2, 1:2, 1:5};

all_perms = perms(1:4);   % 24 × 4

output_dir = fullfile('..', '..', 'bilevel_optim_results', 'job_ioc_results', 'csv_weights');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

for p = 1:size(all_perms, 1)
    perm = all_perms(p, :);

    % Build the row ordering for this permutation by iterating factors
    % in the order given by perm (outermost to innermost)
    row_order = zeros(n_rows, 1);
    k = 0;
    for v1 = factor_ranges{perm(1)}
        for v2 = factor_ranges{perm(2)}
            for v3 = factor_ranges{perm(3)}
                for v4 = factor_ranges{perm(4)}
                    k = k + 1;
                    % Map permuted values back to canonical (subj, leg, phase, speed)
                    vals = zeros(1, 4);
                    vals(perm(1)) = v1;
                    vals(perm(2)) = v2;
                    vals(perm(3)) = v3;
                    vals(perm(4)) = v4;
                    row_order(k) = (vals(1)-1)*20 + (vals(2)-1)*10 + (vals(3)-1)*5 + vals(4);
                end
            end
        end
    end

    % Assemble 41 × 17 cell array (header + 40 data rows)
    C = cell(n_rows + 1, 17);
    C(1, :) = col_headers;
    for r = 1:n_rows
        row = rows(row_order(r));
        C{r+1, 1} = row.label;
        for w = 1:15
            C{r+1, w+1} = row.alpha(w);
        end
        C{r+1, 17} = row.rmse;
    end

    % File name encodes the nesting order (outermost factor first)
    fname = fullfile(output_dir, ...
        sprintf('ioc_weights_%s.csv', strjoin(factor_names(perm), '-')));
    writecell(C, fname);
    fprintf('Written: %s\n', fname);
end

fprintf('\nDone. %d CSV files written to:\n  %s\n', size(all_perms, 1), output_dir);
