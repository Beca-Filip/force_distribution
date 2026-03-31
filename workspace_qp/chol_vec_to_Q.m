function [Q, L] = chol_vec_to_Q(q, n)
%CHOL_VEC_TO_Q Reconstruct a PSD matrix Q from its Cholesky parameterisation.
%
%   [Q, L] = CHOL_VEC_TO_Q(q, n)
%
%   q is a vector of length n*(n+1)/2 containing the lower-triangular entries
%   of the Cholesky factor L, stored in column-major order:
%
%       q = [ L(1,1), L(2,1), ..., L(n,1), L(2,2), L(3,2), ..., L(n,n) ]
%
%   Q is recovered as Q = L * L', which is symmetric PSD by construction.
%
%   For a unique parameterisation, diagonal entries of L (and therefore the
%   corresponding entries in q) should be non-negative.  This is not enforced
%   here but is imposed as a lower bound during IO search (see
%   QP_IO_fmincon_search).
%
%   Inputs:
%     q  - (n*(n+1)/2 x 1) lower-triangular entries of L
%     n  - matrix dimension
%
%   Outputs:
%     Q  - (n x n) symmetric positive-semidefinite matrix
%     L  - (n x n) lower-triangular Cholesky factor

L = zeros(n, n);
L(tril(true(n))) = q;
Q = L * L';

end
