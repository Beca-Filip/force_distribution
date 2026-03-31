function F_invcov = compute_F_invcov(f_data, f_mean)
%COMPUTE_F_INVCOV Compute the whitening (inverse square-root covariance) matrix.
%
%   F_invcov = COMPUTE_F_INVCOV(f_data, f_mean)
%
%   f_data is the force array [n, 1, nsamples, ntrials, nspeeds, nlegs].
%   f_mean is the (n x 1) mean vector from COMPUTE_F_MEAN.
%
%   F_invcov is the (n x n) lower-triangular whitening matrix W such that
%
%       W * Sigma * W' = I
%
%   where Sigma is the sample covariance of the centred forces.  It is
%   obtained as W = inv(L), where L = chol(Sigma, 'lower').  A small Tikhonov
%   regularisation is added to Sigma before factorisation to guard against
%   near-singularity.
%
%   The normalised force used in the QP cost is then:
%       f_norm = F_invcov * (f - f_mean)
%
%   Inputs:
%     f_data  - [n x 1 x nsamples x ntrials x nspeeds x nlegs] force array
%     f_mean  - [n x 1] mean force vector (from compute_f_mean)
%
%   Output:
%     F_invcov - [n x n] whitening transform

n = size(f_data, 1);
F = reshape(f_data, n, []);             % [n x N]
F_centred = F - f_mean;                 % broadcast over columns
N = size(F, 2);

Sigma = F_centred * F_centred' / (N - 1);   % sample covariance [n x n]
Sigma = Sigma + 1e-8 * eye(n);              % Tikhonov regularisation

L = chol(Sigma, 'lower');               % Sigma = L * L'
F_invcov = inv(L);                      % W = L^{-1}  =>  W*Sigma*W' = I

end
