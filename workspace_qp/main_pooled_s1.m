close all; clear all; clc;

% =========================================================================
% Pooled QP inverse optimal control -- SUBJECT 1 (patient 4)
%
% One (Q, l) fitted jointly across every speed and every leg of this subject,
% as opposed to main.m which fits ten separate ones.  This is the
% within-subject half of the question "can a single QP work on a wide set of
% conditions?"; main_pooled_s2.m is the other subject and main_pooled.m
% pools across both.
%
% Results go to main_pooled_s1/:
%   fit.mat                        the pooled fmincon record
%   fit-weights.png, fit-eig.png   the identified Q
%   subject-4-speed-P-leg-L.{mat,png}
%                                  the pooled (Q,l) scored on each condition
%   summary.mat, log.txt, checkpoint.mat
%
% The per-condition .mat files use main.m's variable names, so
%     summarize_main_results('main_pooled_s1')
% audits this folder the same way it audits main/, and the per-condition vs
% pooled comparison (E2) is a join on subject/speed/leg.
%
% RUNTIME.  A pooled objective evaluation solves 10,100 QPs against 1,010 for
% a per-condition fit, so expect roughly ten times the cost of one entry of
% the B5 batch at equal iteration counts.  max_iterations below is the knob
% that bounds it; the run checkpoints every iteration and can be resumed with
%     ck = load('main_pooled_s1/checkpoint.mat');  cfg.theta0 = ck.theta;
% =========================================================================

cfg = struct();
cfg.subject_id  = 4;
cfg.data_path   = fullfile('..', 'Optimization Model Data', 'Patient4.mat');
cfg.results_dir = fullfile(fileparts(mfilename('fullpath')), 'main_pooled_s1');

% All speeds, both legs -- this is what "pooled" means here.
cfg.speed_list  = 1:5;
cfg.leg_list    = [1, 2];
cfg.sample_list = 1:101;
cfg.trial_list  = 1:10;

% Same conditioning as the per-condition fits, so the two are comparable.
cfg.cond_max = 1e4;

% fmincon's default of 1000 iterations is not affordable at this cost.  Raise
% it if the run stops on flag 0 (iteration cap) with a still-decreasing
% objective; the checkpoint lets you resume rather than restart.
cfg.max_iterations = 300;

% l0 is random; fixing the seed makes a multi-hour run repeatable.
cfg.seed = 0;

summary_s1 = run_pooled_fit(cfg);
