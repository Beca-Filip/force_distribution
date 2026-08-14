classdef TestWhitening < matlab.unittest.TestCase
    %TESTWHITENING  Whitening correctness and permutation equivariance.
    %
    %   Run via:  cd workspace_qp; run_qp_tests('whitening')
    %
    %   The point of these tests is one property: with symmetric whitening,
    %   "coordinate i" means "muscle i" and nothing else.  That is what makes
    %   REMAP_QP_TO_SUBJECT exact when carrying a fitted (Q, l) from one
    %   subject to another, since the two subjects list their muscles
    %   differently and share only 33 of them.
    %
    %   test_cholesky_is_not_permutation_equivariant exists to prove the
    %   equivariance tests have teeth: it pins down the failure mode we are
    %   switching away from, so a silent revert to Cholesky cannot pass.
    %
    %   No CasADi needed -- compute_F_invcov and remap_qp_to_subject are
    %   plain MATLAB.

    properties
        f_data      % [n x 1 x ns x nt x nsp x nlg]
        f_mean      % [n x 1]
        n = 8
        perm        % a fixed non-identity permutation
        names       % muscle names in the original order
    end

    methods (TestClassSetup)
        function setup_data(tc)
            rng(20260814);

            ns = 12; nt = 3; nsp = 2; nlg = 2;

            % Correlated forces, so Sigma is genuinely non-diagonal: whitening
            % that mixes coordinates is the only case where order matters.
            N = ns * nt * nsp * nlg;
            mixing = randn(tc.n, tc.n);
            raw    = 5 + mixing * randn(tc.n, N);

            tc.f_data = reshape(raw, tc.n, 1, ns, nt, nsp, nlg);
            tc.f_mean = mean(raw, 2);

            tc.perm  = [5 1 8 3 7 2 6 4]';
            tc.names = arrayfun(@(k) sprintf('m%02d', k), (1:tc.n)', ...
                                'UniformOutput', false);
        end
    end

    methods (Test)

        % ---- basic correctness ------------------------------------------

        function test_symmetric_whitens(tc)
            W = compute_F_invcov(tc.f_data, tc.f_mean, 'symmetric');
            S = tc.sample_cov();
            tc.verifyEqual(W * S * W', eye(tc.n), 'AbsTol', 1e-8, ...
                'W*Sigma*W'' must be the identity.');
        end

        function test_cholesky_whitens(tc)
            W = compute_F_invcov(tc.f_data, tc.f_mean, 'cholesky');
            S = tc.sample_cov();
            tc.verifyEqual(W * S * W', eye(tc.n), 'AbsTol', 1e-8, ...
                'W*Sigma*W'' must be the identity.');
        end

        function test_default_method_is_symmetric(tc)
            W_default = compute_F_invcov(tc.f_data, tc.f_mean);
            W_sym     = compute_F_invcov(tc.f_data, tc.f_mean, 'symmetric');
            tc.verifyEqual(W_default, W_sym, 'AbsTol', 0, ...
                'The default must be symmetric whitening.');
        end

        function test_symmetric_output_is_symmetric_pd(tc)
            W = compute_F_invcov(tc.f_data, tc.f_mean, 'symmetric');
            tc.verifyEqual(W, W', 'AbsTol', 0, 'W must be exactly symmetric.');
            tc.verifyGreaterThan(min(eig(W)), 0, 'W must be positive definite.');
        end

        function test_cholesky_output_is_lower_triangular(tc)
            W = compute_F_invcov(tc.f_data, tc.f_mean, 'cholesky');
            tc.verifyEqual(norm(triu(W, 1), 'fro'), 0, 'AbsTol', 1e-12, ...
                'Cholesky whitening must stay lower triangular.');
        end

        function test_unknown_method_errors(tc)
            tc.verifyError(@() compute_F_invcov(tc.f_data, tc.f_mean, 'nope'), ...
                'compute_F_invcov:unknownMethod');
        end

        % ---- the property that matters -----------------------------------

        function test_symmetric_is_permutation_equivariant(tc)
            % W(P*Sigma*P') must equal P*W(Sigma)*P'.
            p = tc.perm;
            W_orig = compute_F_invcov(tc.f_data, tc.f_mean, 'symmetric');
            W_perm = compute_F_invcov(tc.f_data(p, :, :, :, :, :), ...
                                      tc.f_mean(p), 'symmetric');
            tc.verifyEqual(W_perm, W_orig(p, p), 'AbsTol', 1e-9, ...
                ['Symmetric whitening must commute with relabelling the ' ...
                 'muscles; otherwise coordinate i is not muscle i.']);
        end

        function test_cholesky_is_not_permutation_equivariant(tc)
            % Guard test: pins the failure mode we are moving away from, so a
            % revert to Cholesky cannot slip past the test above.
            p = tc.perm;
            W_orig = compute_F_invcov(tc.f_data, tc.f_mean, 'cholesky');
            W_perm = compute_F_invcov(tc.f_data(p, :, :, :, :, :), ...
                                      tc.f_mean(p), 'cholesky');
            discrepancy = max(abs(W_perm(:) - reshape(W_orig(p, p), [], 1)));
            tc.verifyGreaterThan(discrepancy, 1e-6, ...
                ['Cholesky whitening is expected to be order dependent. If ' ...
                 'this passes, the fixture stopped being correlated and the ' ...
                 'equivariance test above proves nothing.']);
        end

        function test_normalised_forces_follow_the_muscle(tc)
            % The user-facing consequence: f_norm coordinate i belongs to
            % muscle i, so relabelling just permutes f_norm.
            p = tc.perm;
            W_a = compute_F_invcov(tc.f_data, tc.f_mean, 'symmetric');
            W_b = compute_F_invcov(tc.f_data(p, :, :, :, :, :), ...
                                   tc.f_mean(p), 'symmetric');

            f = tc.f_data(:, 1, 4, 2, 1, 2);
            fn_a = W_a * (f      - tc.f_mean);
            fn_b = W_b * (f(p)   - tc.f_mean(p));

            tc.verifyEqual(fn_b, fn_a(p), 'AbsTol', 1e-9, ...
                'Normalised forces must permute with the muscle labels.');
        end

        % ---- end to end: whitening + remap --------------------------------

        function test_remapped_cost_matches_across_orderings(tc)
            % The claim the whole switch is for: fit (Q, l) in subject A's
            % ordering, carry it to subject B's ordering with
            % remap_qp_to_subject, and the QP cost of the SAME physical force
            % vector is unchanged.  This only holds because the whitening is
            % equivariant -- with Cholesky the two f_norm differ and so does
            % the cost.
            p = tc.perm;
            names_a = tc.names;
            names_b = tc.names(p);

            W_a = compute_F_invcov(tc.f_data, tc.f_mean, 'symmetric');
            W_b = compute_F_invcov(tc.f_data(p, :, :, :, :, :), ...
                                   tc.f_mean(p), 'symmetric');

            % An arbitrary PSD Q and linear term in A's ordering.
            L    = tril(randn(tc.n) + 2 * eye(tc.n));
            Q_a  = L * L';
            l_a  = randn(tc.n, 1);

            [Q_b, l_b] = remap_qp_to_subject(Q_a, l_a, names_a, names_b);

            f    = tc.f_data(:, 1, 7, 1, 2, 1);
            fn_a = W_a * (f    - tc.f_mean);
            fn_b = W_b * (f(p) - tc.f_mean(p));

            cost_a = 0.5 * fn_a' * Q_a * fn_a + l_a' * fn_a;
            cost_b = 0.5 * fn_b' * Q_b * fn_b + l_b' * fn_b;

            tc.verifyEqual(cost_b, cost_a, 'RelTol', 1e-9, ...
                'A remapped QP must score the same force identically.');
        end

        function test_remap_is_exact_on_a_pure_relabelling(tc)
            % With no muscles added or dropped the remap must be the exact
            % permutation -- no fill values involved.
            p = tc.perm;
            L   = tril(randn(tc.n) + 2 * eye(tc.n));
            Q_a = L * L';
            l_a = randn(tc.n, 1);

            [Q_b, l_b, info] = remap_qp_to_subject(Q_a, l_a, tc.names, tc.names(p));

            tc.verifyEqual(Q_b, Q_a(p, p), 'AbsTol', 0);
            tc.verifyEqual(l_b, l_a(p),    'AbsTol', 0);
            tc.verifyEqual(info.n_inserted, 0);
            tc.verifyEqual(info.n_dropped,  0);
        end
    end

    methods (Access = private)
        function S = sample_cov(tc)
            % Must mirror compute_F_invcov, regularisation included.
            F  = reshape(tc.f_data, tc.n, []);
            Fc = F - tc.f_mean;
            S  = Fc * Fc' / (size(F, 2) - 1);
            S  = (S + S') / 2 + 1e-8 * eye(tc.n);
        end
    end
end
