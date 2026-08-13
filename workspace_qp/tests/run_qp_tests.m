function results = run_qp_tests(varargin)
%RUN_QP_TESTS  Entry point for the workspace_qp test suites.
%
%   results = RUN_QP_TESTS()               runs every suite
%   results = RUN_QP_TESTS('gradient')     runs TestMultiConditionGradient
%   results = RUN_QP_TESTS('synergy')      runs TestSynergyPipeline
%
%   Run from the workspace_qp directory:
%     cd workspace_qp
%     results = run_qp_tests();
%     disp(results)
%
%   WHY THIS SCRIPT EXISTS
%   ----------------------
%   Two path traps make it unsafe to call runtests directly on these files.
%
%   1. MATLAB binds function names inside a classdef when the CLASS IS FIRST
%      LOADED, not when its methods run.  A test class that calls addpath in
%      its own TestClassSetup is therefore too late: any function it resolved
%      at load time keeps the old binding.  This bit hard with rmse -- the
%      project's utils/error_utils/rmse flattens with a(:) and returns a
%      scalar, while MATLAB's built-in rmse reduces along dim 1 and returns
%      an ARRAY.  With the built-in bound, QP_IO_inner_loop silently returns
%      a non-scalar objective.  Paths must be set BEFORE the class loads.
%
%   2. runtests('tests/TestX') fails to build a suite from a relative path in
%      -batch mode; matlab.unittest.TestSuite.fromFile works reliably.
%
%   Both are handled here, so this is the supported way to run the tests.

% ---- Paths (set before any test class is loaded) ------------------------
casadi_path = 'C:/Users/filip/Documents/GitHub/casadi-3.7.0-windows64-matlab2018b';

this_dir      = fileparts(mfilename('fullpath'));   % workspace_qp/tests
workspace_dir = fileparts(this_dir);                % workspace_qp
repo_dir      = fileparts(workspace_dir);           % repo root

addpath(casadi_path);
addpath(workspace_dir);
addpath(fullfile(repo_dir, 'workspace_do'));            % eq/ineq_constraint_function
addpath(fullfile(repo_dir, 'utils', 'error_utils'));    % scalar rmse

% ---- Guard: the project rmse must win over MATLAB's built-in ------------
% Checked here, at script scope, where name resolution is re-evaluated.
if ~isequal(size(rmse(ones(2,2,2), zeros(2,2,2))), [1 1])
    error(['run_qp_tests:wrongRmse  MATLAB''s built-in rmse is shadowing ' ...
           'the project version in %s. The objective would be an array ' ...
           'instead of a scalar.'], fullfile(repo_dir, 'utils', 'error_utils'));
end

% ---- Select suites ------------------------------------------------------
if isempty(varargin)
    which_suite = 'all';
else
    which_suite = lower(varargin{1});
end

% fromFile returns matlab.unittest.Test elements, so the accumulator must be
% Test.empty -- TestSuite.empty cannot be concatenated with them.
suite = matlab.unittest.Test.empty(1, 0);

if any(strcmp(which_suite, {'all', 'gradient'}))
    suite = [suite, matlab.unittest.TestSuite.fromFile( ...
        fullfile(this_dir, 'TestMultiConditionGradient.m'))];
end

if any(strcmp(which_suite, {'all', 'synergy'}))
    suite = [suite, matlab.unittest.TestSuite.fromFile( ...
        fullfile(this_dir, 'TestSynergyPipeline.m'))];
end

if isempty(suite)
    error('run_qp_tests:unknownSuite', 'Unknown suite selector "%s".', which_suite);
end

% ---- Run ----------------------------------------------------------------
results = run(suite);

fprintf('\n==== SUMMARY ====\n');
fprintf('passed=%d failed=%d incomplete=%d  (%.1f s)\n', ...
    sum([results.Passed]), sum([results.Failed]), ...
    sum([results.Incomplete]), sum([results.Duration]));

end
