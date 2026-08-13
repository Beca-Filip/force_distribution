classdef TestMultiConditionGradient < matlab.unittest.TestCase
%TESTMULTICONDITIONGRADIENT  Verifies that QP_IO_inner_loop returns a correct
%                            objective AND gradient when several speeds and
%                            legs are pooled into a single fit.
%
%   Background
%   ----------
%   QP_IO_inner_loop originally computed its gradient as
%
%       M  = ns * nt;
%       r  = reshape(Fout(:,1,:,:,1,1) - Fref(:,1,:,:,1,1), n, M);
%       dF = reshape(dFout(:,:,:,:,1,1), n, ntheta, M);
%
%   which silently discards every condition except (speed_list(1), leg_list(1)).
%   The objective E was already pooled correctly (rmse flattens with a(:)),
%   so the mismatch is invisible unless the gradient is checked directly with
%   nspeeds > 1 or nlegs > 1.  The pre-existing finite-difference test in
%   TestSynergyPipeline only ever ran with speed_list = 1, leg_list = 1, so it
%   could not catch this.  Every test below is run with 2 speeds x 2 legs.
%
%   Run from MATLAB with:
%     cd workspace_qp
%     results = tests/run_qp_tests('gradient');
%
%   Use run_qp_tests -- do NOT call runtests on this file directly.  MATLAB
%   binds function names inside a classdef when the class is first loaded, so
%   any addpath performed by this class's own setup is already too late (see
%   run_qp_tests.m for the full explanation).

    properties (Constant)
        CASADI_PATH = 'C:/Users/filip/Documents/GitHub/casadi-3.7.0-windows64-matlab2018b'
    end

    properties
        data
        model
        vars
        n
        theta0
        sample_list
        trial_list
        speed_list
        leg_list
    end

    % =====================================================================
    methods (TestClassSetup)

        % Setup is driven from a SINGLE TestClassSetup method, because MATLAB
        % does not guarantee the execution order of multiple TestClassSetup
        % methods (splitting these made the data builder run before addpath).
        %
        % The three steps must live in SEPARATE method scopes: MATLAB binds
        % function names per function scope, so calls made in the same
        % function that ran addpath still resolve against the OLD path.  With
        % everything inlined here, rmse bound to MATLAB's built-in (which
        % reduces along dim 1) instead of the project's flattening version.
        function setup_environment(tc)
            tc.add_required_paths();
            tc.assert_project_rmse_is_active();
            tc.build_multi_condition_data();
        end

    end

    % =====================================================================
    methods (Access = private)

        function add_required_paths(tc)
            % Safety net for a session where run_qp_tests was not used.
            % Paths are derived from this file's location, not pwd: the test
            % framework runs the class with tests/ as the working folder.
            this_dir      = fileparts(mfilename('fullpath'));  % workspace_qp/tests
            workspace_dir = fileparts(this_dir);               % workspace_qp
            repo_dir      = fileparts(workspace_dir);          % repo root

            addpath(tc.CASADI_PATH);
            addpath(workspace_dir);
            addpath(fullfile(repo_dir, 'workspace_do'));          % eq/ineq_constraint_function
            addpath(fullfile(repo_dir, 'utils', 'error_utils'));  % scalar rmse
        end

        function assert_project_rmse_is_active(tc)
            % rmse must be the project's flattening version, not MATLAB's
            % built-in.  The built-in reduces along dim 1 and would make
            % QP_IO_inner_loop return an ARRAY objective instead of a scalar.
            tc.assertEqual(rmse(ones(2,2,2), zeros(2,2,2)), 1, 'AbsTol', 1e-12, ...
                ['rmse resolved to MATLAB''s built-in instead of the project ' ...
                 'utility. Launch the tests via tests/run_qp_tests.m, which ' ...
                 'sets the path before this class is loaded.']);
        end

        function build_multi_condition_data(tc)
            % Synthetic data with genuinely different content in every
            % (speed, leg) cell, so that each condition makes a distinct
            % contribution to both the objective and the gradient.
            %
            % n = 10 against the ne = 5 equality constraints hardcoded in
            % form_casadi_qp_model leaves 5 free directions, so the QP
            % solution is not pinned down by the constraints alone and the
            % cost weights genuinely influence it.
            n_muscles = 10;
            n_samples = 4;
            n_trials  = 2;
            n_speeds  = 2;
            n_legs    = 2;
            ne        = 5;

            rng(7);

            f = abs(randn(n_muscles, 1, n_samples, n_trials, n_speeds, n_legs)) + 0.5;

            % Loose bounds so solutions sit strictly inside the box and the
            % active set is just the equality constraints.
            tc.data.f    = f;
            tc.data.fmin = 0.2 * f;
            tc.data.fmax = 3.0 * f;

            % Distinct A per (speed, leg); b = A*f keeps every QP feasible.
            a_data = zeros(ne, n_muscles, n_samples, n_trials, n_speeds, n_legs);
            b_data = zeros(ne, 1,         n_samples, n_trials, n_speeds, n_legs);
            for sp = 1:n_speeds
                for lg = 1:n_legs
                    a_cond = randn(ne, n_muscles);   % differs per speed/leg
                    for kk = 1:n_samples
                        for tt = 1:n_trials
                            a_data(:, :, kk, tt, sp, lg) = a_cond;
                            b_data(:, :, kk, tt, sp, lg) = a_cond * f(:, :, kk, tt, sp, lg);
                        end
                    end
                end
            end
            tc.data.A = a_data;
            tc.data.b = b_data;

            tc.data.f_mean   = compute_f_mean(f);
            tc.data.F_invcov = compute_F_invcov(f, tc.data.f_mean);

            [tc.model, tc.vars] = form_casadi_qp_model(n_muscles);

            % Tight solver tolerances: OSQP's loose defaults would otherwise
            % dominate the finite-difference error and make the FD checks
            % below unreliable.
            sol_opt = struct;
            sol_opt.osqp.verbose  = 0;
            sol_opt.osqp.eps_abs  = 1e-10;
            sol_opt.osqp.eps_rel  = 1e-10;
            sol_opt.osqp.max_iter = 200000;
            tc.model.solver('osqp', sol_opt);

            tc.n           = n_muscles;
            tc.sample_list = 1:n_samples;
            tc.trial_list  = 1:n_trials;
            tc.speed_list  = 1:n_speeds;
            tc.leg_list    = 1:n_legs;

            % theta0: L = I (so Q = I), l random but fixed
            nq = n_muscles * (n_muscles + 1) / 2;
            mask = tril(true(n_muscles));
            [rr, cc] = find(mask);
            q0 = zeros(nq, 1);
            q0(rr == cc) = 1;
            tc.theta0 = [q0; 0.3 * randn(n_muscles, 1)];
        end

    end

    % =====================================================================
    % Shape / sanity
    % =====================================================================
    methods (Test)

        function test_gradient_is_finite_and_correct_length(tc)
            [E, dE] = tc.eval_pooled(tc.theta0);
            tc.verifyTrue(isscalar(E) && isfinite(E) && E >= 0, ...
                'E must be a finite non-negative scalar.');
            tc.verifyEqual(numel(dE), numel(tc.theta0));
            tc.verifyTrue(all(isfinite(dE)), 'dE must be finite.');
        end

        function test_gradient_is_nonzero(tc)
            % A zero gradient would make the FD comparisons vacuous.
            [~, dE] = tc.eval_pooled(tc.theta0);
            tc.verifyGreaterThan(norm(dE), 1e-8, ...
                'Gradient is numerically zero -- the test problem is degenerate.');
        end

    end

    % =====================================================================
    % Finite-difference checks with multiple speeds AND legs
    % =====================================================================
    methods (Test)

        function test_gradient_fd_multi_speed_multi_leg(tc)
            % The central test: this fails against the pre-fix implementation,
            % which built its gradient from (speed 1, leg 1) only.
            theta = tc.theta0;
            [~, dE_an] = tc.eval_pooled(theta);

            rng(3);
            nq = tc.n * (tc.n + 1) / 2;
            % Sample from both halves of theta: Cholesky entries and l entries.
            check_idx = [randperm(nq, 5), nq + randperm(tc.n, 3)];

            for ci = check_idx
                dE_fd = tc.central_difference(theta, ci);
                tc.verify_gradient_component(dE_an(ci), dE_fd, ci);
            end
        end

        function test_gradient_fd_single_condition_regression(tc)
            % The pre-fix code path was correct for one speed/leg; confirm the
            % fix did not disturb it.
            theta = tc.theta0;
            [~, dE_an] = QP_IO_inner_loop(theta, tc.data, tc.vars, tc.model, ...
                tc.sample_list, tc.trial_list, 1, 1);

            rng(4);
            nq = tc.n * (tc.n + 1) / 2;
            check_idx = [randperm(nq, 3), nq + randperm(tc.n, 2)];

            for ci = check_idx
                dE_fd = tc.central_difference(theta, ci, 1, 1);
                tc.verify_gradient_component(dE_an(ci), dE_fd, ci);
            end
        end

    end

    % =====================================================================
    % Exact decomposition identity (no finite differences involved)
    % =====================================================================
    methods (Test)

        function test_pooled_objective_matches_condition_decomposition(tc)
            % E = sqrt(sum_c S_c / sum_c N_c),  S_c = E_c^2 * N_c
            [E_pool, ~] = tc.eval_pooled(tc.theta0);

            s_total = 0;
            n_total = 0;
            for sp = tc.speed_list
                for lg = tc.leg_list
                    [E_c, N_c] = tc.eval_condition(tc.theta0, sp, lg);
                    s_total = s_total + E_c^2 * N_c;
                    n_total = n_total + N_c;
                end
            end
            E_expected = sqrt(s_total / n_total);

            tc.verifyEqual(E_pool, E_expected, 'RelTol', 1e-9, ...
                'Pooled objective does not decompose over conditions.');
        end

        function test_pooled_gradient_matches_condition_decomposition(tc)
            % Each QP is solved independently, so the pooled residual sum is
            % the sum of the per-condition residual sums:
            %
            %   g_c      = E_c * N_c * dE_c          (unnormalised sum r.*dF)
            %   dE_pool  = (1/(E*N)) * sum_c g_c
            %
            % This is an exact algebraic identity -- it pins down the pooled
            % gradient to solver precision without any finite differencing,
            % and it fails loudly if conditions are dropped or misaligned.
            [E_pool, dE_pool] = tc.eval_pooled(tc.theta0);

            g_total = zeros(numel(tc.theta0), 1);
            n_total = 0;
            for sp = tc.speed_list
                for lg = tc.leg_list
                    [E_c, N_c, dE_c] = tc.eval_condition(tc.theta0, sp, lg);
                    g_total = g_total + E_c * N_c * dE_c;
                    n_total = n_total + N_c;
                end
            end
            dE_expected = g_total / (E_pool * n_total);

            tc.verifyEqual(dE_pool, dE_expected, 'RelTol', 1e-7, 'AbsTol', 1e-10, ...
                'Pooled gradient does not decompose over conditions.');
        end

        function test_gradient_invariant_to_condition_ordering(tc)
            % Both E and dE are sums over conditions, so reversing the order
            % of speed_list / leg_list must not change either.
            %
            % This is a sharp test for the dropped-condition bug and needs no
            % finite differencing: code that reads only list POSITION 1 picks
            % a different physical condition once the list is reversed, so its
            % answer moves.  Correct code is exactly invariant.
            [E_fwd, dE_fwd] = QP_IO_inner_loop(tc.theta0, tc.data, tc.vars, ...
                tc.model, tc.sample_list, tc.trial_list, ...
                tc.speed_list, tc.leg_list);

            [E_rev, dE_rev] = QP_IO_inner_loop(tc.theta0, tc.data, tc.vars, ...
                tc.model, tc.sample_list, tc.trial_list, ...
                fliplr(tc.speed_list), fliplr(tc.leg_list));

            tc.verifyEqual(E_rev, E_fwd, 'RelTol', 1e-9, ...
                'Pooled objective depends on the order of speed_list/leg_list.');
            tc.verifyEqual(dE_rev, dE_fwd, 'RelTol', 1e-7, 'AbsTol', 1e-10, ...
                ['Pooled gradient depends on the order of speed_list/leg_list ' ...
                 '-- conditions are being dropped or misaligned.']);
        end

    end

    % =====================================================================
    % Helpers
    % =====================================================================
    methods (Access = private)

        function [E, dE] = eval_pooled(tc, theta)
            [E, dE] = QP_IO_inner_loop(theta, tc.data, tc.vars, tc.model, ...
                tc.sample_list, tc.trial_list, tc.speed_list, tc.leg_list);
        end

        function [E_c, N_c, dE_c] = eval_condition(tc, theta, sp, lg)
            [E_c, dE_c] = QP_IO_inner_loop(theta, tc.data, tc.vars, tc.model, ...
                tc.sample_list, tc.trial_list, sp, lg);
            N_c = tc.n * numel(tc.sample_list) * numel(tc.trial_list);
        end

        function dE_fd = central_difference(tc, theta, idx, sp, lg)
            if nargin < 4
                sp = tc.speed_list;
                lg = tc.leg_list;
            end
            h = 1e-6;
            theta_plus       = theta;  theta_plus(idx)  = theta_plus(idx)  + h;
            theta_minus      = theta;  theta_minus(idx) = theta_minus(idx) - h;

            E_plus  = QP_IO_inner_loop(theta_plus,  tc.data, tc.vars, tc.model, ...
                tc.sample_list, tc.trial_list, sp, lg);
            E_minus = QP_IO_inner_loop(theta_minus, tc.data, tc.vars, tc.model, ...
                tc.sample_list, tc.trial_list, sp, lg);

            dE_fd = (E_plus - E_minus) / (2 * h);
        end

        function verify_gradient_component(tc, dE_an, dE_fd, idx)
            % Relative where the component is large, absolute where it is
            % near zero (FD noise floor is set by the QP solver tolerance).
            tol = max(1e-6, 1e-3 * abs(dE_fd));
            tc.verifyEqual(dE_an, dE_fd, 'AbsTol', tol, ...
                sprintf(['Gradient mismatch at theta(%d): analytic = %.9f, ' ...
                         'finite difference = %.9f (tol %.2e)'], ...
                        idx, dE_an, dE_fd, tol));
        end

    end

end
