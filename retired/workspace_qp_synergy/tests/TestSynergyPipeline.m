classdef TestSynergyPipeline < matlab.unittest.TestCase
%TESTSYNERGYPIPELINE  Unit and integration tests for the synergy-constrained
%                     QP inverse-optimal-control pipeline.
%
%   Run from MATLAB with:
%     cd workspace_qp
%     results = runtests('tests/TestSynergyPipeline');
%     disp(results)

    properties (Constant)
        CASADI_PATH = 'C:/Users/filip/Documents/GitHub/casadi-3.7.0-windows64-matlab2018b'
    end

    properties
        muscle_names_p4
        muscle_names_p5
        synergy_groups
        syn_info_p4
        syn_info_p5
        small_data         % minimal data struct for fast integration tests
        small_model
        small_vars
        small_syn_info
    end

    % =====================================================================
    methods (TestClassSetup)

        function add_paths(tc)
            % Paths are derived from this file's location rather than pwd:
            % the test framework runs the class with tests/ as the working
            % folder, so pwd-relative paths resolve one level too deep.
            repo = tc.repo_root();
            addpath(tc.CASADI_PATH);
            addpath(fullfile(repo, 'workspace_qp'));            % shared QP helpers
            addpath(fullfile(repo, 'retired', 'workspace_qp_synergy'));  % retired synergy sources
            addpath(fullfile(repo, 'workspace_do'));            % eq_constraint_function
            addpath(fullfile(repo, 'utils', 'error_utils'));    % scalar rmse
        end

        function load_muscle_names(tc)
            repo = tc.repo_root();
            mn4 = load(fullfile(repo, 'patient_4_muscle_names.mat'));
            mn5 = load(fullfile(repo, 'patient_5_muscle_names.mat'));
            tc.muscle_names_p4 = mn4.patient_4_muscle_names;
            tc.muscle_names_p5 = mn5.patient_5_muscle_names;
            tc.synergy_groups  = { ...
                {'glmax1','glmax2','glmax3','vasmed','vasint','vaslat'}, ...
                {'soleus','gasmed','gaslat'}, ...
                {'iliacus','psoas','recfem'}, ...
                {'semimem','semiten','bflh','bfsh'}, ...
                {'tibant','edl'}, ...
                {'glmed1','glmed2','glmed3','glmin1','glmin2','glmin3', ...
                 'addbrev','addlong','addmagDist','addmagIsch','addmagMid','addmagProx'}, ...
            };
            tc.syn_info_p4 = build_synergy_info(tc.muscle_names_p4, tc.synergy_groups);
            tc.syn_info_p5 = build_synergy_info(tc.muscle_names_p5, tc.synergy_groups);
        end

        function build_small_data(tc)
            % Build a minimal synthetic data struct for fast integration tests.
            % n=6 muscles, 2 synergy groups (size 4 and 2), 1 singleton.
            n          = 6;
            nsamples   = 5;
            ntrials    = 2;
            nspeeds    = 1;
            nlegs      = 1;

            rng(0);
            f = abs(randn(n, 1, nsamples, ntrials, nspeeds, nlegs)) + 0.1;
            tc.small_data.f    = f;
            tc.small_data.fmin = 0.8 * f;
            tc.small_data.fmax = 1.2 * f;

            ne = 5;  % must match ne hardcoded in form_casadi_qp_model
            A  = randn(ne, n);
            % Set b = A*f per sample/trial so the true force always satisfies
            % the equality constraint, guaranteeing QP feasibility.
            b_data = zeros(ne, 1, nsamples, ntrials, nspeeds, nlegs);
            for kk = 1:nsamples
                for tt = 1:ntrials
                    b_data(:, :, kk, tt, 1, 1) = A * f(:, :, kk, tt, 1, 1);
                end
            end
            tc.small_data.A = repmat(A, [1 1 nsamples ntrials nspeeds nlegs]);
            tc.small_data.b = b_data;

            tc.small_data.f_mean   = compute_f_mean(f);
            tc.small_data.F_invcov = compute_F_invcov(f, tc.small_data.f_mean);

            small_muscle_names = {'g1','g1b','g1c','g1d','g2','g2b'};
            small_synergy_groups = {{'g1'},{'g2'}};
            tc.small_syn_info = build_synergy_info(small_muscle_names, small_synergy_groups);

            [tc.small_model, tc.small_vars] = form_casadi_qp_model(n);
            sol_opt = struct; sol_opt.osqp.verbose = 0;
            tc.small_model.solver('osqp', sol_opt);
        end

    end

    % =====================================================================
    % Tests: build_synergy_info
    % =====================================================================
    methods (Test)

        function test_group_count_p4(tc)
            si = tc.syn_info_p4;
            tc.verifyEqual(si.n, numel(tc.muscle_names_p4));
        end

        function test_group_count_p5(tc)
            si = tc.syn_info_p5;
            tc.verifyEqual(si.n, numel(tc.muscle_names_p5));
        end

        function test_no_muscle_in_two_groups(tc)
            % Each muscle index may appear in at most one group
            all_idx = [];
            for k = 1:tc.syn_info_p4.n_groups
                all_idx = [all_idx, tc.syn_info_p4.group_indices{k}]; %#ok<AGROW>
            end
            tc.verifyEqual(numel(all_idx), numel(unique(all_idx)), ...
                'A muscle appears in more than one synergy group.');
        end

        function test_full_q_indices_unique(tc)
            si = tc.syn_info_p4;
            tc.verifyEqual(numel(si.full_q_indices), numel(unique(si.full_q_indices)), ...
                'Duplicate entries in full_q_indices.');
        end

        function test_full_q_indices_in_range(tc)
            si = tc.syn_info_p4;
            tc.verifyTrue(all(si.full_q_indices >= 1) && all(si.full_q_indices <= si.n_full_q));
        end

        function test_n_theta_q_formula(tc)
            % n_theta_q should equal sum of m_k*(m_k+1)/2 + n_singletons
            si = tc.syn_info_p4;
            expected = 0;
            for k = 1:si.n_groups
                m = numel(si.group_indices{k});
                expected = expected + m*(m+1)/2;
            end
            expected = expected + numel(si.singleton_indices);
            tc.verifyEqual(si.n_theta_q, expected);
        end

        function test_diag_mask_count(tc)
            % Number of diagonal entries = sum of group sizes + n_singletons = n
            si = tc.syn_info_p4;
            tc.verifyEqual(sum(si.diag_mask), si.n);
        end

        function test_p5_recfem_in_exactly_one_group(tc)
            % recfem should appear in exactly one group for P5
            si = tc.syn_info_p5;
            names = tc.muscle_names_p5;
            recfem_idx = find(strcmpi(names, 'recfem'));
            tc.verifyFalse(isempty(recfem_idx), 'recfem not found in P5 muscle names.');
            count = 0;
            for k = 1:si.n_groups
                if any(si.group_indices{k} == recfem_idx)
                    count = count + 1;
                end
            end
            tc.verifyEqual(count, 1, 'recfem should belong to exactly one synergy group.');
        end

    end

    % =====================================================================
    % Tests: synergy_theta_to_Q
    % =====================================================================
    methods (Test)

        function test_Q_is_symmetric(tc)
            si = tc.syn_info_p4;
            theta = tc.make_theta0(si);
            [Q, ~] = synergy_theta_to_Q(theta, si);
            tc.verifyTrue(issymmetric(Q), 'Q must be symmetric.');
        end

        function test_Q_is_psd(tc)
            si = tc.syn_info_p4;
            theta = tc.make_theta0(si);
            [Q, ~] = synergy_theta_to_Q(theta, si);
            lam = eig((Q+Q')/2);
            tc.verifyTrue(all(lam >= -1e-10), 'Q must be PSD.');
        end

        function test_Q_sparsity_pattern(tc)
            % Off-synergy Q entries must be exactly zero
            si = tc.syn_info_p4;
            theta = tc.make_theta0(si);
            [Q, ~] = synergy_theta_to_Q(theta, si);

            % Build expected sparsity mask
            n = si.n;
            allowed = false(n, n);
            for k = 1:si.n_groups
                gi = si.group_indices{k};
                allowed(gi, gi) = true;
            end
            for s = si.singleton_indices
                allowed(s, s) = true;
            end

            off_block_max = max(abs(Q(~allowed)));
            tc.verifyTrue(isempty(off_block_max) || off_block_max == 0, ...
                sprintf('Non-zero off-synergy Q entry: %.2e', off_block_max));
        end

        function test_roundtrip_l(tc)
            si = tc.syn_info_p4;
            l_in = randn(si.n, 1);
            q0   = zeros(si.n_theta_q, 1); q0(si.diag_mask) = 1;
            [~, ~, ~, l_out] = synergy_theta_to_Q([q0; l_in], si);
            tc.verifyEqual(l_out, l_in, 'AbsTol', 1e-14);
        end

    end

    % =====================================================================
    % Tests: QP_synergy_IO_inner_loop
    % =====================================================================
    methods (Test)

        function test_inner_loop_finite_E(tc)
            [E, dE] = tc.run_inner_loop();
            tc.verifyTrue(isfinite(E) && E >= 0, 'E must be finite and non-negative.');
            tc.verifyTrue(all(isfinite(dE)), 'dE must be finite.');
        end

        function test_inner_loop_gradient_length(tc)
            [~, dE] = tc.run_inner_loop();
            tc.verifyEqual(numel(dE), tc.small_syn_info.n_theta);
        end

        function test_gradient_finite_difference_check(tc)
            % Verify analytical gradient against central finite differences
            % on a 3-element perturbation of theta (cheap).
            si    = tc.small_syn_info;
            theta = tc.make_theta0(si);
            h     = 1e-5;

            [~, dE_an] = QP_synergy_IO_inner_loop(theta, si, ...
                tc.small_data, tc.small_vars, tc.small_model, 1:5, 1:2, 1, 1);

            % Check 3 random parameter directions
            rng(1);
            check_idx = randperm(si.n_theta, min(3, si.n_theta));
            for ci = check_idx
                tp = theta; tp(ci) = tp(ci) + h;
                tm = theta; tm(ci) = tm(ci) - h;
                Ep = QP_synergy_IO_inner_loop(tp, si, ...
                    tc.small_data, tc.small_vars, tc.small_model, 1:5, 1:2, 1, 1);
                Em = QP_synergy_IO_inner_loop(tm, si, ...
                    tc.small_data, tc.small_vars, tc.small_model, 1:5, 1:2, 1, 1);
                dE_fd = (Ep - Em) / (2*h);
                tc.verifyEqual(dE_an(ci), dE_fd, 'AbsTol', 1e-4, ...
                    sprintf('Gradient mismatch at theta(%d): analytic=%.6f  FD=%.6f', ...
                        ci, dE_an(ci), dE_fd));
            end
        end

    end

    % =====================================================================
    % Tests: fmincon search (smoke test — just checks it runs and converges)
    % =====================================================================
    methods (Test)

        function test_fmincon_smoke(tc)
            si    = tc.small_syn_info;
            theta0 = tc.make_theta0(si);

            % Use very loose tolerances and few iterations for speed
            % (we only test that it runs without error and returns finite output)
            orig_path = path;
            cleanup   = onCleanup(@() path(orig_path));

            % Override fmincon options via a wrapper that tightens limits
            fun = @(th) QP_synergy_IO_inner_loop(th, si, ...
                tc.small_data, tc.small_vars, tc.small_model, 1:5, 1:2, 1, 1);

            lb = -inf(si.n_theta, 1); lb(si.diag_mask) = 0;
            opts = optimoptions(@fmincon, ...
                'SpecifyObjectiveGradient', true, ...
                'Display', 'off', ...
                'MaxIterations', 10, ...
                'MaxFunctionEvaluations', 50);

            [th_opt, fval] = fmincon(fun, theta0, [], [], [], [], lb, [], [], opts);

            tc.verifyTrue(isfinite(fval), 'fmincon returned non-finite objective.');
            tc.verifyTrue(all(isfinite(th_opt)), 'fmincon returned non-finite theta.');
        end

    end

    % =====================================================================
    % Helpers
    % =====================================================================
    methods (Access = private)

        function repo = repo_root(~)
            % <repo>/retired/workspace_qp_synergy/tests/this_file.m  ->  <repo>
            this_dir = fileparts(mfilename('fullpath'));
            repo     = fileparts(fileparts(fileparts(this_dir)));
        end

        function theta = make_theta0(~, si)
            q0 = zeros(si.n_theta_q, 1); q0(si.diag_mask) = 1;
            rng(42);
            theta = [q0; randn(si.n, 1)];
        end

        function [E, dE] = run_inner_loop(tc)
            si    = tc.small_syn_info;
            theta = tc.make_theta0(si);
            [E, dE] = QP_synergy_IO_inner_loop(theta, si, ...
                tc.small_data, tc.small_vars, tc.small_model, 1:5, 1:2, 1, 1);
        end

    end

end
