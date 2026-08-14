function F_invcov = compute_F_invcov(f_data, f_mean, varargin)
%COMPUTE_F_INVCOV Compute the whitening (inverse square-root covariance) matrix.
%
%   F_invcov = COMPUTE_F_INVCOV(f_data, f_mean)
%   F_invcov = COMPUTE_F_INVCOV(f_data, f_mean, method)
%
%   f_data is the force array [n, 1, nsamples, ntrials, nspeeds, nlegs].
%   f_mean is the (n x 1) mean vector from COMPUTE_F_MEAN.
%
%   F_invcov is the (n x n) whitening matrix W such that
%
%       W * Sigma * W' = I
%
%   where Sigma is the sample covariance of the centred forces.  The
%   normalised force used in the QP cost is then
%
%       f_norm = F_invcov * (f - f_mean)
%
%   method is 'symmetric' (default) or 'cholesky'.
%
%
%   WHY SYMMETRIC IS THE DEFAULT
%   ----------------------------
%   A whitening matrix is only defined up to a left orthogonal factor: if
%   W*Sigma*W' = I then so does (R*W)*Sigma*(R*W)' for any orthogonal R.  The
%   choice of R decides what "coordinate i" means, and that choice matters as
%   soon as a fitted Q is moved between subjects.
%
%   'cholesky' takes W = inv(chol(Sigma,'lower')), which is lower triangular.
%   Coordinate i is then muscle i's residual after regressing out muscles
%   1..i-1, so its meaning depends on the ORDER the muscles happen to be
%   listed in.  Formally, Cholesky whitening is not permutation-equivariant:
%   for a permutation P,
%
%       W_chol(P*Sigma*P')  ~=  P * W_chol(Sigma) * P'.
%
%   The two subjects here make that concrete.  Patient 4 has 35 muscles,
%   patient 5 has 34, sharing 33; tfl exists only in patient 5 and sits at
%   its index 18.  Under Cholesky whitening, inserting tfl there changes the
%   conditioning set of EVERY later coordinate, so patient 4's coordinate 25
%   and patient 5's coordinate 25 are not the same quantity even though they
%   are the same muscle.  A Q fitted on one subject then cannot be carried to
%   the other by re-indexing, no matter how careful the re-indexing is --
%   which is precisely what REMAP_QP_TO_SUBJECT tries to do.
%
%   'symmetric' takes the symmetric inverse square root W = Sigma^(-1/2),
%   computed from the eigendecomposition.  Every matrix function of a
%   symmetric matrix is permutation-equivariant,
%
%       Sigma^(-1/2) evaluated at P*Sigma*P'  =  P * Sigma^(-1/2) * P',
%
%   because P*Sigma*P' = (P*V)*D*(P*V)' is an eigendecomposition with the same
%   eigenvalues D.  Coordinate i therefore refers to muscle i regardless of
%   listing order, and REMAP_QP_TO_SUBJECT becomes exact on the shared block
%   rather than approximate.  (Sigma^(-1/2) is also the whitening matrix
%   closest to the identity in Frobenius norm, so it mixes muscles as little
%   as any whitening can.)
%
%   Switching the whitening changes the coordinates the cost is expressed in,
%   so it INVALIDATES any Q, l fitted under the previous convention.  Refits
%   are required.  'cholesky' is retained only to reproduce those older runs.
%
%   Inputs:
%     f_data  - [n x 1 x nsamples x ntrials x nspeeds x nlegs] force array
%     f_mean  - [n x 1] mean force vector (from compute_f_mean)
%     method  - 'symmetric' (default) | 'cholesky'
%
%   Output:
%     F_invcov - [n x n] whitening transform; symmetric positive definite for
%                'symmetric', lower triangular for 'cholesky'

if nargin < 3 || isempty(varargin{1})
    method = 'symmetric';
else
    method = lower(varargin{1});
end

reg = 1e-8;                             % Tikhonov regularisation of Sigma

n = size(f_data, 1);
F = reshape(f_data, n, []);             % [n x N]
F_centred = F - f_mean;                 % broadcast over columns
N = size(F, 2);

Sigma = F_centred * F_centred' / (N - 1);   % sample covariance [n x n]
Sigma = (Sigma + Sigma') / 2;               % kill round-off asymmetry
Sigma = Sigma + reg * eye(n);               % Tikhonov regularisation

switch method
    case 'symmetric'
        % W = Sigma^(-1/2) = V * diag(d.^-0.5) * V'.
        % eig on a symmetric input returns an orthonormal V, and the result
        % is independent of eigenvalue ordering or of the basis chosen inside
        % a repeated eigenvalue's eigenspace, so W is unique.
        [V, D] = eig(Sigma, 'vector');
        D = max(D, reg);                % guard: clip round-off negatives
        F_invcov = V * diag(1 ./ sqrt(D)) * V';
        F_invcov = (F_invcov + F_invcov') / 2;   % exact symmetry

    case 'cholesky'
        L = chol(Sigma, 'lower');       % Sigma = L * L'
        F_invcov = inv(L);              % W = L^{-1}  =>  W*Sigma*W' = I

    otherwise
        error('compute_F_invcov:unknownMethod', ...
              'method must be ''symmetric'' or ''cholesky'', got "%s".', method);
end

end
