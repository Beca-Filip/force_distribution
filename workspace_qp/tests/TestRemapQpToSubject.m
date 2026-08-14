classdef TestRemapQpToSubject < matlab.unittest.TestCase
%TESTREMAPQPTOSUBJECT  Tests for remap_qp_to_subject.
%
%   Covers the three structural cases (matched / dropped / inserted) against
%   independently computed references, plus the algebraic properties the
%   remap must satisfy: PSD preservation, the eigenvalue union, exactness of
%   permutation and round-trip, and homogeneity under rescaling of (Q, l).
%
%   The real patient 4 / patient 5 name lists are used where the point is
%   that the function handles the actual 35-vs-34 mismatch.
%
%   Run from MATLAB with:
%     cd workspace_qp
%     results = tests/run_qp_tests('remap');
%
%   Use run_qp_tests -- see run_qp_tests.m for why calling runtests directly
%   on this file is unreliable.

    properties
        names_p4
        names_p5
        Q_p4          % a well-conditioned PSD matrix on the P4 muscle set
        l_p4
    end

    % =====================================================================
    methods (TestClassSetup)

        function setup_environment(tc)
            tc.add_required_paths();
            tc.load_real_muscle_names();
        end

    end

    methods (Access = private)

        function add_required_paths(~)
            this_dir      = fileparts(mfilename('fullpath'));  % workspace_qp/tests
            workspace_dir = fileparts(this_dir);               % workspace_qp
            addpath(workspace_dir);
        end

        function load_real_muscle_names(tc)
            this_dir = fileparts(mfilename('fullpath'));
            repo_dir = fileparts(fileparts(this_dir));

            m4 = load(fullfile(repo_dir, 'patient_4_muscle_names.mat'));
            m5 = load(fullfile(repo_dir, 'patient_5_muscle_names.mat'));
            f4 = fieldnames(m4);
            f5 = fieldnames(m5);
            tc.names_p4 = cellstr(m4.(f4{1}));
            tc.names_p5 = cellstr(m5.(f5{1}));
            tc.names_p4 = tc.names_p4(:);
            tc.names_p5 = tc.names_p5(:);

            % A PSD, comfortably conditioned Q on the P4 set.
            n = numel(tc.names_p4);
            rng(11);
            A = randn(n, n);
            tc.Q_p4 = A * A' / n + eye(n);     % strictly positive definite
            tc.l_p4 = randn(n, 1);
        end

        function [Q, l] = random_psd(~, n, seed)
            rng(seed);
            A = randn(n, n);
            Q = A * A' / n + 0.5 * eye(n);
            l = randn(n, 1);
        end

    end

    % =====================================================================
    % Structural cases
    % =====================================================================
    methods (Test)

        function test_identity_remap_is_exact(tc)
            % Same names in the same order must be a no-op, bit for bit.
            [Q_to, l_to, info] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p4);

            tc.verifyEqual(Q_to, tc.Q_p4, 'AbsTol', 0, ...
                'Identity remap altered Q.');
            tc.verifyEqual(l_to, tc.l_p4, 'AbsTol', 0, ...
                'Identity remap altered l.');
            tc.verifyEqual(info.n_matched, numel(tc.names_p4));
            tc.verifyEqual(info.n_dropped, 0);
            tc.verifyEqual(info.n_inserted, 0);
        end

        function test_pure_permutation_matches_explicit_permutation(tc)
            % names_to = names_from(perm)  =>  Q_to = Q_from(perm, perm).
            n = numel(tc.names_p4);
            rng(21);
            perm = randperm(n)';

            names_shuffled = tc.names_p4(perm);
            [Q_to, l_to, info] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, names_shuffled);

            tc.verifyEqual(info.n_matched, n);
            tc.verifyEqual(info.n_dropped, 0);
            tc.verifyEqual(info.n_inserted, 0);
            tc.verifyEqual(Q_to, tc.Q_p4(perm, perm), 'AbsTol', 1e-12, ...
                'Permuted remap does not equal the explicit permutation of Q.');
            tc.verifyEqual(l_to, tc.l_p4(perm), 'AbsTol', 1e-12, ...
                'Permuted remap does not equal the explicit permutation of l.');
        end

        function test_drop_only_gives_principal_submatrix(tc)
            % Removing names must yield exactly the principal submatrix.
            keep_idx = [1:10, 12:numel(tc.names_p4)]';   % drop muscle 11
            names_kept = tc.names_p4(keep_idx);

            [Q_to, l_to, info] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, names_kept);

            tc.verifyEqual(info.n_dropped, 1);
            tc.verifyEqual(info.n_inserted, 0);
            tc.verifyEqual(info.names_dropped{1}, tc.names_p4{11});
            tc.verifyEqual(Q_to, tc.Q_p4(keep_idx, keep_idx), 'AbsTol', 1e-12);
            tc.verifyEqual(l_to, tc.l_p4(keep_idx), 'AbsTol', 1e-12);
        end

        function test_inserted_muscle_is_diagonal_and_uncoupled(tc)
            % An inserted muscle gets a positive diagonal and zero coupling.
            names_extended = [tc.names_p4; {'brand_new_muscle'}];

            [Q_to, l_to, info] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, names_extended, ...
                'q_fill', 'median', 'l_fill', 'zero');

            tc.verifyEqual(info.n_inserted, 1);
            u = info.idx_to_inserted;

            off_diag_row = Q_to(u, :);
            off_diag_row(u) = [];
            tc.verifyEqual(norm(off_diag_row), 0, 'AbsTol', 0, ...
                'Inserted muscle must have zero coupling to all others.');

            tc.verifyEqual(Q_to(u, u), median(diag(tc.Q_p4)), 'AbsTol', 1e-12, ...
                'Inserted diagonal should be the median of the matched diagonal.');
            tc.verifyEqual(l_to(u), 0, 'AbsTol', 0, ...
                'Default l_fill is zero.');

            % The original block must be untouched.
            tc.verifyEqual(Q_to(1:end-1, 1:end-1), tc.Q_p4, 'AbsTol', 1e-12);
        end

    end

    % =====================================================================
    % Algebraic properties
    % =====================================================================
    methods (Test)

        function test_eigenvalues_are_union_of_block_and_fills(tc)
            % Q_to is similar to blkdiag(Q_matched, diag(fills)), so its
            % spectrum is exactly the union of the two.  This is the precise
            % statement of "the fill creates no new null direction".
            [Q_to, ~, info] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5, 'q_fill', 'median');

            Q_matched = tc.Q_p4(info.idx_from_matched, info.idx_from_matched);
            expected  = sort([eig(Q_matched); info.q_fill_values]);
            actual    = sort(eig(Q_to));

            tc.verifyEqual(actual, expected, 'RelTol', 1e-9, 'AbsTol', 1e-10, ...
                'Spectrum of the remapped Q is not the expected union.');
        end

        function test_psd_is_preserved_on_real_name_sets(tc)
            [Q_to, ~, info] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5);

            tc.verifyTrue(info.is_psd, 'Remapped Q reported as not PSD.');
            tc.verifyGreaterThan(min(eig(Q_to)), 0, ...
                'Remapped Q should stay positive definite here.');
        end

        function test_scale_homogeneity(tc)
            % The QP minimiser is invariant to (Q,l) -> (cQ, cl), so the
            % remap must commute with rescaling -- otherwise the fill would
            % silently change the relative weighting.
            c = 7.5;
            [Q_a, l_a] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5, 'l_fill', 'median');
            [Q_b, l_b] = remap_qp_to_subject( ...
                c * tc.Q_p4, c * tc.l_p4, tc.names_p4, tc.names_p5, 'l_fill', 'median');

            tc.verifyEqual(Q_b, c * Q_a, 'RelTol', 1e-12, ...
                'Remap does not commute with rescaling of Q.');
            tc.verifyEqual(l_b, c * l_a, 'RelTol', 1e-12, ...
                'Remap does not commute with rescaling of l.');
        end

        function test_symmetry_is_enforced(tc)
            [Q_to, ~] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5);
            tc.verifyEqual(Q_to, Q_to', 'AbsTol', 0, ...
                'Remapped Q must be exactly symmetric.');
        end

    end

    % =====================================================================
    % Real patient 4 <-> patient 5 behaviour
    % =====================================================================
    methods (Test)

        function test_real_p4_to_p5_membership(tc)
            [~, ~, info] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5);

            tc.verifyEqual(info.n_matched,  33);
            tc.verifyEqual(info.n_dropped,   2);
            tc.verifyEqual(info.n_inserted,  1);
            tc.verifyEqual(sort(lower(info.names_dropped)),  {'edl'; 'fdl'});
            tc.verifyEqual(lower(info.names_inserted), {'tfl'});
        end

        function test_real_p5_to_p4_membership(tc)
            n5 = numel(tc.names_p5);
            [Q5, l5] = tc.random_psd(n5, 31);

            [~, ~, info] = remap_qp_to_subject(Q5, l5, tc.names_p5, tc.names_p4);

            tc.verifyEqual(info.n_matched,  33);
            tc.verifyEqual(info.n_dropped,   1);
            tc.verifyEqual(info.n_inserted,  2);
            tc.verifyEqual(lower(info.names_dropped), {'tfl'});
            tc.verifyEqual(sort(lower(info.names_inserted)), {'edl'; 'fdl'});
        end

        function test_round_trip_preserves_the_common_block(tc)
            % P4 -> P5 -> P4.  The 33 shared muscles must come back bit for
            % bit; edl/fdl cannot (they were dropped) and are refilled.
            [Q5, l5] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5);
            [Q4b, l4b, info_back] = remap_qp_to_subject( ...
                Q5, l5, tc.names_p5, tc.names_p4);

            common = intersect(lower(tc.names_p4), lower(tc.names_p5));
            [~, idx_orig] = ismember(common, lower(tc.names_p4));
            [~, idx_back] = ismember(common, lower(tc.names_p4));

            tc.verifyEqual(Q4b(idx_back, idx_back), tc.Q_p4(idx_orig, idx_orig), ...
                'AbsTol', 1e-12, ...
                'Round trip corrupted the shared 33x33 block of Q.');
            tc.verifyEqual(l4b(idx_back), tc.l_p4(idx_orig), 'AbsTol', 1e-12, ...
                'Round trip corrupted the shared entries of l.');

            % edl and fdl were genuinely lost and must have been re-inserted.
            tc.verifyEqual(sort(lower(info_back.names_inserted)), {'edl'; 'fdl'});
        end

        function test_donor_fill_uses_named_donor_diagonal(tc)
            % tfl (hip abductor) filled from glmed1 rather than a statistic.
            [Q_to, l_to, info] = remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5, ...
                'q_fill', 'donor', 'l_fill', 'donor', ...
                'donor_map', {'tfl', 'glmed1'});

            donor_idx = find(strcmpi(tc.names_p4, 'glmed1'), 1);
            u = info.idx_to_inserted;

            tc.verifyEqual(Q_to(u, u), tc.Q_p4(donor_idx, donor_idx), 'AbsTol', 1e-12);
            tc.verifyEqual(l_to(u), tc.l_p4(donor_idx), 'AbsTol', 1e-12);

            % Still uncoupled: donor supplies the diagonal only, never a row.
            off_diag_row = Q_to(u, :);
            off_diag_row(u) = [];
            tc.verifyEqual(norm(off_diag_row), 0, 'AbsTol', 0, ...
                'Donor fill must not copy the donor row (that risks singularity).');
        end

    end

    % =====================================================================
    % Input validation
    % =====================================================================
    methods (Test)

        function test_duplicate_names_error(tc)
            bad = tc.names_p4;
            bad{2} = bad{1};
            tc.verifyError(@() remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, bad, tc.names_p5), ...
                'remap_qp_to_subject:duplicateNames');
        end

        function test_size_mismatch_error(tc)
            tc.verifyError(@() remap_qp_to_subject( ...
                tc.Q_p4(1:10, 1:10), tc.l_p4, tc.names_p4, tc.names_p5), ...
                'remap_qp_to_subject:sizeMismatch');
        end

        function test_no_overlap_error(tc)
            foreign = arrayfun(@(k) sprintf('zzz_%d', k), 1:5, 'UniformOutput', false);
            tc.verifyError(@() remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, foreign), ...
                'remap_qp_to_subject:noOverlap');
        end

        function test_non_positive_q_fill_error(tc)
            tc.verifyError(@() remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5, 'q_fill', 0), ...
                'remap_qp_to_subject:nonPositiveFill');
        end

        function test_missing_donor_map_error(tc)
            tc.verifyError(@() remap_qp_to_subject( ...
                tc.Q_p4, tc.l_p4, tc.names_p4, tc.names_p5, 'q_fill', 'donor'), ...
                'remap_qp_to_subject:missingDonorMap');
        end

    end

end
