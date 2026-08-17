close all; clear all; clc;

% =========================================================================
% Pooled QP inverse optimal control -- SUBJECT 2 (patient 5)
%
% Identical to main_pooled_s1.m in every respect except the subject: one
% (Q, l) fitted jointly across every speed and every leg.  Keeping the two
% scripts identical apart from subject_id / data_path / results_dir is
% deliberate -- any difference between the two fits is then a difference
% between the subjects, not between the runs.
%
% Note this subject has 34 muscles against subject 4's 35 (P4 has edl and
% fdl, P5 has tfl), so the two Q matrices live in different spaces and are
% not directly comparable entry by entry.  remap_qp_to_subject.m is what
% bridges them, and main_pooled.m is where that happens.
%
% Results go to main_pooled_s2/, same layout as main_pooled_s1/.  Audit with
%     summarize_main_results('main_pooled_s2')
%
% See main_pooled_s1.m for the runtime and resume notes.
% =========================================================================

cfg = struct();
cfg.subject_id  = 5;
cfg.data_path   = fullfile('..', 'Optimization Model Data', 'Patient5.mat');
cfg.results_dir = fullfile(fileparts(mfilename('fullpath')), 'main_pooled_s2');

cfg.speed_list  = 1:5;
cfg.leg_list    = [1, 2];
cfg.sample_list = 1:101;
cfg.trial_list  = 1:10;

cfg.cond_max = 1e4;

cfg.max_iterations = 300;

cfg.seed = 0;

summary_s2 = run_pooled_fit(cfg);
