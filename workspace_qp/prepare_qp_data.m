function data = prepare_qp_data(data)
%PREPARE_QP_DATA  Repair the force bounds and add the normalisation fields.
%
%   data = PREPARE_QP_DATA(data)
%
%   Two jobs, both of which must be done IDENTICALLY for the per-condition
%   fits (main.m) and the pooled ones (main_pooled_*.m).  If they diverge the
%   per-condition vs pooled comparison compares data preparation as much as it
%   compares objectives, so this lives in one place and both call it.
%
%   1. BOUND REPAIR.  The QP needs the observed force to be strictly feasible:
%      fmin < f < fmax.  Where the recorded bounds touch or cross the observed
%      force the QP would be infeasible, or feasible only at a vertex, so the
%      offending bound is pushed out by 0.3%.
%
%   2. NORMALISATION.  f_mean and F_invcov are computed from ALL of the
%      subject's data (every sample, trial, speed and leg), not from the
%      subset being fitted.  That is deliberate: the normalisation defines the
%      coordinate system in which Q lives, and it must not change between a
%      per-condition fit and a pooled one, or the two Q matrices would not be
%      expressed in the same coordinates and could not be compared.
%
%      F_invcov defaults to the symmetric whitening Sigma^(-1/2), which is
%      permutation-equivariant and therefore safe to remap between subjects.
%      See COMPUTE_F_INVCOV.
%
%   See also COMPUTE_F_MEAN, COMPUTE_F_INVCOV, MAIN, RUN_POOLED_FIT.

% ---- 1. Bound repair ----------------------------------------------------
data.fmax(data.fmax <= data.f) = 1.003 * data.f(data.fmax <= data.f);
data.fmin(data.fmin >= data.f) = 0.997 * data.f(data.fmin >= data.f);

% The two lines below are inherited verbatim from the original script and are
% inert on the present data: after the two assignments above, every one of the
% three masks is empty for both patients.  They are kept for exact fidelity
% with the fits already produced, NOT because they are right -- the first one
% indexes the left side with |fmax - f| < eps and the right side with
% |fmax - fmin| < eps, so if it ever did fire on a data set where those two
% masks select different numbers of elements, it would error.  Fix it if a new
% subject makes it fire; do not silently change it now, because that would
% change every fit.
epsil = 0.01;
data.fmax(abs(data.fmax - data.f) < epsil) = ...
    data.fmax(abs(data.fmax - data.fmin) < epsil) + epsil;
data.fmin(abs(data.f - data.fmin) < epsil) = max( ...
    zeros(size(data.fmin(abs(data.f - data.fmin) < epsil))), ...
    data.fmin(abs(data.f - data.fmin) < epsil) - epsil);

% ---- 2. Normalisation ---------------------------------------------------
data.f_mean   = compute_f_mean(data.f);
data.F_invcov = compute_F_invcov(data.f, data.f_mean);

end
