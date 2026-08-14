classdef TestConditioning < matlab.unittest.TestCase
    %TESTCONDITIONING  The cond(Q) <= cond_max guarantee and the trace gauge fix.
    %
    %   Run via:  cd workspace_qp; run_qp_tests('conditioning')
    %
    %   Covers THETA_TO_QL and QP_TRACE_CONSTRAINT. The claim being pinned is
    %
    %       Q = L*L' + eps*I,  eps = n/cond_max,  trace(Q) = n
    %       =>  lambda_min >= eps,  lambda_max <= n,  so cond(Q) <= cond_max
    %
    %   and, separately, that the trace EQUALITY is what removes the
    %   (Q,l) -> (cQ, cl) scale gauge.
    %
    %   Two tests exist to give the others teeth:
    %     test_gauge_freedom_is_real          -- shows the gauge exists at all,
    %                                            by solving the DOC twice
    %     test_bound_fails_without_constraint -- shows cond blows past cond_max
    %                                            once the sphere is left
    %   Without those, the guarantee tests could pass on a vacuous fixture.
    %
    %   No CasADi needed except where noted; the DOC gauge demo uses quadprog.

    properties
        n = 12
        cond_max = 1e4
        nq
        ntheta
        eps_shift
    end

    methods (TestClassSetup)
        function setup_sizes(tc)
            rng(20260814);
            tc.nq        = tc.n * (tc.n + 1) / 2;
            tc.ntheta    = tc.nq + tc.n;
            tc.eps_shift = tc.n / tc.cond_max;
        end
    end

    methods (Test)

        % ---- theta_to_Ql ------------------------------------------------

        function test_shift_is_n_over_cond_max(tc)
            [~, ~, ~, eps_shift] = theta_to_Ql(tc.random_theta(), tc.n, tc.cond_max);
            tc.verifyEqual(eps_shift, tc.n / tc.cond_max, 'RelTol', 1e-15);
        end

        function test_Q_equals_LLt_plus_shift(tc)
            theta = tc.random_theta();
            [Q, ~, L, eps_shift] = theta_to_Ql(theta, tc.n, tc.cond_max);
            tc.verifyEqual(Q, L*L' + eps_shift*eye(tc.n), 'AbsTol', 1e-12, ...
                'Q must be exactly L*L'' plus the shift.');
        end

        function test_L_is_unshifted_and_lower_triangular(tc)
            % kkt_sensitivity uses L for dQ/dL, which must NOT see the shift.
            theta = tc.random_theta();
            [~, ~, L] = theta_to_Ql(theta, tc.n, tc.cond_max);
            tc.verifyEqual(norm(triu(L, 1), 'fro'), 0, 'AbsTol', 1e-14);
            L_ref = chol_vec_to_Q_factor(theta(1:tc.nq), tc.n);
            tc.verifyEqual(L, L_ref, 'AbsTol', 0, ...
                'L must be the raw Cholesky factor, shift-free.');
        end

        function test_lambda_min_at_least_shift(tc)
            % True for ANY theta, including a rank-deficient L.
            q = zeros(tc.nq, 1);                 % L = 0  =>  Q = eps*I
            theta = [q; randn(tc.n, 1)];
            [Q, ~, ~, eps_shift] = theta_to_Ql(theta, tc.n, tc.cond_max);
            tc.verifyGreaterThanOrEqual(min(eig(Q)), eps_shift - 1e-12, ...
                'The shift must floor the spectrum even when L is singular.');
        end

        function test_Q_is_exactly_symmetric(tc)
            [Q, ~] = theta_to_Ql(tc.random_theta(), tc.n, tc.cond_max);
            tc.verifyEqual(Q, Q', 'AbsTol', 0);
        end

        function test_default_cond_max_is_1e4(tc)
            theta = tc.random_theta();
            [Q_default, ~, ~, eps_default] = theta_to_Ql(theta, tc.n);
            [Q_explicit, ~, ~, eps_explicit] = theta_to_Ql(theta, tc.n, 1e4);
            tc.verifyEqual(Q_default,   Q_explicit,   'AbsTol', 0);
            tc.verifyEqual(eps_default, eps_explicit, 'AbsTol', 0);
        end

        function test_inf_cond_max_disables_shift(tc)
            theta = tc.random_theta();
            [Q, ~, L, eps_shift] = theta_to_Ql(theta, tc.n, Inf);
            tc.verifyEqual(eps_shift, 0);
            tc.verifyEqual(Q, L*L', 'AbsTol', 1e-12, ...
                'cond_max = Inf must reproduce the legacy unshifted Q.');
        end

        function test_l_block_is_untouched(tc)
            theta = tc.random_theta();
            [~, l] = theta_to_Ql(theta, tc.n, tc.cond_max);
            tc.verifyEqual(l, reshape(theta(tc.nq+1:end), tc.n, 1), 'AbsTol', 0);
        end

        function test_bad_theta_length_errors(tc)
            tc.verifyError(@() theta_to_Ql(zeros(tc.ntheta+1, 1), tc.n, tc.cond_max), ...
                'theta_to_Ql:badTheta');
        end

        function test_bad_cond_max_errors(tc)
            tc.verifyError(@() theta_to_Ql(tc.random_theta(), tc.n, -1), ...
                'theta_to_Ql:badCondMax');
            tc.verifyError(@() theta_to_Ql(tc.random_theta(), tc.n, [1 2]), ...
                'theta_to_Ql:badCondMax');
        end

        % ---- qp_trace_constraint -----------------------------------------

        function test_residual_matches_trace_directly(tc)
            theta = tc.random_theta();
            [~, ceq] = qp_trace_constraint(theta, tc.n, tc.cond_max);
            Q = theta_to_Ql(theta, tc.n, tc.cond_max);
            tc.verifyEqual(ceq, trace(Q) - tc.n, 'AbsTol', 1e-10, ...
                'The residual must equal trace(Q) - n.');
        end

        function test_feasible_init_is_exactly_feasible(tc)
            % L0 = sqrt(1-eps)*I, the initialisation main.m uses.
            theta0 = tc.feasible_theta();
            [~, ceq] = qp_trace_constraint(theta0, tc.n, tc.cond_max);
            tc.verifyEqual(ceq, 0, 'AbsTol', 1e-12, ...
                'main.m''s starting point must satisfy trace(Q) = n exactly.');

            Q0 = theta_to_Ql(theta0, tc.n, tc.cond_max);
            tc.verifyEqual(Q0, eye(tc.n), 'AbsTol', 1e-12, ...
                'That starting point should give exactly Q0 = I.');
        end

        function test_constraint_gradient_matches_finite_differences(tc)
            theta = tc.random_theta();
            [~, ~, ~, gceq] = qp_trace_constraint(theta, tc.n, tc.cond_max);

            h = 1e-6;
            gceq_fd = zeros(tc.ntheta, 1);
            for k = 1:tc.ntheta
                tp = theta; tp(k) = tp(k) + h;
                tm = theta; tm(k) = tm(k) - h;
                [~, ceq_p] = qp_trace_constraint(tp, tc.n, tc.cond_max);
                [~, ceq_m] = qp_trace_constraint(tm, tc.n, tc.cond_max);
                gceq_fd(k) = (ceq_p - ceq_m) / (2*h);
            end
            tc.verifyEqual(gceq, gceq_fd, 'AbsTol', 1e-6);
        end

        function test_constraint_gradient_shape_and_l_block(tc)
            % fmincon wants [numel(theta) x n_ceq]; l must not be constrained,
            % since an unconstrained l is what lets the trace fix the gauge.
            [c, ceq, gc, gceq] = qp_trace_constraint(tc.random_theta(), tc.n, tc.cond_max);
            tc.verifyEmpty(c);
            tc.verifyEmpty(gc);
            tc.verifySize(ceq,  [1 1]);
            tc.verifySize(gceq, [tc.ntheta 1]);
            tc.verifyEqual(gceq(tc.nq+1:end), zeros(tc.n, 1), 'AbsTol', 0);
        end

        % ---- the guarantee ------------------------------------------------

        function test_cond_bound_holds_on_the_constraint_surface(tc)
            % The whole point: anywhere on the sphere, cond(Q) <= cond_max.
            for trial = 1:200
                theta = tc.feasible_theta_random();
                Q = theta_to_Ql(theta, tc.n, tc.cond_max);
                e = eig(Q);
                tc.verifyLessThanOrEqual(max(e)/min(e), tc.cond_max * (1 + 1e-9), ...
                    sprintf('cond bound violated on draw %d.', trial));
            end
        end

        function test_cond_bound_holds_for_a_rank_one_L(tc)
            % Worst case for the old parameterisation: L*L' has rank 1, so the
            % unshifted Q is singular (cond = Inf). The shift must rescue it.
            v = randn(tc.n, 1);
            L = tril(v * [1, zeros(1, tc.n-1)]);   % only column 1 populated
            q = L(tril(true(tc.n)));
            q = q * sqrt(tc.n * (1 - tc.eps_shift)) / norm(q);   % onto the sphere
            theta = [q; randn(tc.n, 1)];

            Q = theta_to_Ql(theta, tc.n, tc.cond_max);
            e = eig(Q);
            tc.verifyLessThanOrEqual(max(e)/min(e), tc.cond_max * (1 + 1e-9));

            Q_unshifted = theta_to_Ql(theta, tc.n, Inf);
            tc.verifyGreaterThan(cond(Q_unshifted), 1e12, ...
                'Fixture check: the unshifted Q here should be near-singular.');
        end

        function test_bound_fails_without_constraint(tc)
            % Teeth: off the sphere, or with the shift disabled, the bound is
            % not just unproven but actually violated. If this ever passes,
            % the guarantee tests above are vacuous.
            v = randn(tc.n, 1);
            L = tril(v * [1, zeros(1, tc.n-1)]);
            q = L(tril(true(tc.n)));
            theta_off = [q * 1e3; randn(tc.n, 1)];   % far outside the sphere

            Q_off = theta_to_Ql(theta_off, tc.n, tc.cond_max);
            tc.verifyGreaterThan(cond(Q_off), tc.cond_max, ...
                ['Leaving the trace surface must break the bound; if it does ' ...
                 'not, the constraint is not doing any work.']);
        end

        % ---- fmincon integration -------------------------------------------

        function test_fmincon_converges_onto_the_trace_surface(tc)
            % Verifies qp_trace_constraint is correctly FORMED for fmincon --
            % residual sign, gradient orientation ([ntheta x 1], not a row),
            % and the l block left free.  A transposed or mis-signed gradient
            % still "runs" but converges to the wrong place or not at all.
            %
            % Deliberately uses a simple smooth objective rather than the QP
            % objective: this is a test of the constraint, and the IO problem
            % itself does not converge in any budget suitable for a unit test
            % (see test_fmincon_maintains_the_cond_bound).
            theta_start = tc.feasible_theta();
            target      = 0.3;

            opts = optimoptions(@fmincon, ...
                'SpecifyObjectiveGradient',  true, ...
                'SpecifyConstraintGradient', true, ...
                'Display',                   'off', ...
                'ConstraintTolerance',       1e-6, ...
                'MaxIterations',             200);

            [theta_end, ~, exitflag] = fmincon( ...
                @(th) quartic_objective(th, target), theta_start, ...
                [], [], [], [], [], [], ...
                @(th) qp_trace_constraint(th, tc.n, tc.cond_max), opts);

            tc.verifyGreaterThan(exitflag, 0, 'fmincon must converge.');

            [~, ceq] = qp_trace_constraint(theta_end, tc.n, tc.cond_max);
            tc.verifyLessThanOrEqual(abs(ceq), 1e-6, ...
                'Converged run must sit on the trace surface.');

            Q_end = theta_to_Ql(theta_end, tc.n, tc.cond_max);
            e     = eig(Q_end);
            tc.verifyEqual(trace(Q_end), tc.n, 'AbsTol', 1e-6);
            tc.verifyLessThanOrEqual(max(e)/min(e), tc.cond_max * (1 + 1e-6), ...
                'The guarantee must hold at the converged point.');

            % The objective pulls every entry toward `target`, which off the
            % surface would give ||q||^2 = nq*target^2. It must not get there.
            tc.verifyNotEqual(theta_end(1:tc.nq), target * ones(tc.nq, 1));
        end

        % ---- the gauge ----------------------------------------------------

        function test_gauge_freedom_is_real(tc)
            % Motivation check, solved directly: rescaling (Q,l) -> (cQ, cl)
            % leaves the DOC minimiser unchanged, so the IO objective has an
            % exact flat direction. This is what the trace equality removes.
            m = 5;
            B = randn(m); Q = B*B' + eye(m);
            l = randn(m, 1);
            A = randn(2, m);
            f_feasible = randn(m, 1);
            b = A * f_feasible;
            lb = f_feasible - 3;
            ub = f_feasible + 3;

            opts = optimoptions('quadprog', 'Display', 'off', ...
                                'OptimalityTolerance', 1e-12);
            f_1 = quadprog(Q,     l,     [], [], A, b, lb, ub, [], opts);
            c   = 7.3;
            f_c = quadprog(c*Q,   c*l,   [], [], A, b, lb, ub, [], opts);

            tc.verifyEqual(f_c, f_1, 'AbsTol', 1e-7, ...
                'The DOC minimiser must be invariant to a positive rescale.');
        end

        function test_trace_constraint_blocks_the_gauge(tc)
            % And the constraint is what stops the optimiser from moving along
            % that flat direction: a rescaled theta is infeasible.
            theta = tc.feasible_theta();
            [~, ceq_1] = qp_trace_constraint(theta, tc.n, tc.cond_max);
            tc.verifyEqual(ceq_1, 0, 'AbsTol', 1e-12);

            for c = [0.5, 2.0, 7.3]
                theta_scaled = theta;
                theta_scaled(1:tc.nq)     = sqrt(c) * theta(1:tc.nq);   % Q -> cQ
                theta_scaled(tc.nq+1:end) = c       * theta(tc.nq+1:end);
                [~, ceq_c] = qp_trace_constraint(theta_scaled, tc.n, tc.cond_max);
                tc.verifyGreaterThan(abs(ceq_c), 1e-3, ...
                    sprintf('Rescaling by c = %g must leave the feasible set.', c));
            end
        end
    end

    methods (Access = private)
        function theta = random_theta(tc)
            theta = randn(tc.ntheta, 1);
        end

        function theta = feasible_theta(tc)
            % main.m's initialisation: L0 = sqrt(1-eps)*I  =>  Q0 = I.
            L0 = sqrt(1 - tc.eps_shift) * eye(tc.n);
            theta = [L0(tril(true(tc.n))); randn(tc.n, 1)];
        end

        function theta = feasible_theta_random(tc)
            % A random point projected onto the sphere ||q||^2 = n(1-eps).
            q = randn(tc.nq, 1);
            q = q * sqrt(tc.n * (1 - tc.eps_shift)) / norm(q);
            theta = [q; randn(tc.n, 1)];
        end
    end
end

% -------------------------------------------------------------------------
function [f, g] = quartic_objective(theta, target)
% Smooth, strictly convex, minimised at theta = target everywhere -- so the
% only thing keeping the iterate off that point is the trace constraint.
f = sum((theta - target).^4);
g = 4 * (theta - target).^3;
end

% -------------------------------------------------------------------------
function L = chol_vec_to_Q_factor(q, n)
% Independent reimplementation of the factor, so the test does not simply
% restate chol_vec_to_Q's own output.
L = zeros(n);
L(tril(true(n))) = q;
end
