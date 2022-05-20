function solio = solve_io_kkt_model_normalized(modelio, varargin)
%SOLVE_IO_KKT_MODEL_NORMALIZED solves the Keshavaraz et al. (2011) inverse
%optimization formulation model as a constrained LSQ problem, using lsqlin.
%
%   modelio = SOLVE_IO_KKT_MODEL_NORMALIZED(modelio)
%   modelio = SOLVE_IO_KKT_MODEL_NORMALIZED(modelio, lsqopt)
%   lsqopt ~ options for lsqlin

% If options are passed
if nargin > 1
    lsqopt = varargin{1};
% Define default options
else
    lsqopt = optimoptions(@lsqlin, ...
                          'Display', 'Iter', ...
                          'ConstraintTolerance', 1e-6, ...
                          'MaxIterations', 500);
end

% Extract lists
sample_list = modelio.sample_list;
trial_list = modelio.trial_list;
speed_list = modelio.speed_list;
leg_list = modelio.leg_list;

% Define matrices
C = modelio.regressor;
d = zeros(size(modelio.regressor, 1), 1);

A = [];
b = [];

Aeq = modelio.normalization_cf_Aeq;
beq = modelio.normalization_cf_beq;

lb = modelio.lower_bound;
ub = [];

x0 = [];

% Solve
[v_sol,resnorm,residual,~,~,~] = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,varargin);

% Memorize vector solution
solio.resnorm = resnorm;
solio.residual = residual;
solio.v_sol = v_sol;

% Extract subvectors
solio.alpha = v_sol(1 : modelio.num_cf);

% Prealocate multipliers
solio.lambda = zeros(modelio.num_ec, 1, modelio.num_samples, modelio.num_trials, modelio.num_speeds, modelio.num_legs);
solio.mu = cell(modelio.num_samples, modelio.num_trials, modelio.num_speeds, modelio.num_legs);

% Extract multipliers from block matrix
% Counters for the loop
cntsamples = 1;
cnttrials = 1;
cntspeeds = 1;
cntlegs = 1;
% Loop over trials
for trial = trial_list
    % Reset speed counter
    cntspeeds = 1;
    % Loop over speeds
    for speed = speed_list
        % Reset leg counter
        cntlegs = 1;
        % Loop over legs
        for leg = leg_list
            % Reset sample counter
            cntsamples = 1;
            % Loop over samples
            for k = sample_list
                
                % Get indices for current equality constraint multipliers
                [~, ic] = indexing_function(modelio, 'ec_multipliers', k, trial, speed, leg);
                % Assign the value to the the matrix
                solio.lambda(:, :, cntsamples, cnttrials, cntspeeds, cntlegs) = v_sol(ic);
                
                % Get indices for current equality constraint multipliers
                [~, ic] = indexing_function(modelio, 'aic_multipliers', k, trial, speed, leg);
                % Assign the value to the the matrix
                solio.mu{cntsamples, cnttrials, cntspeeds, cntlegs} = v_sol(ic);
                                
                % Augment sample counter
                cntsamples = cntsamples + 1;
            end
            % Augment leg count
            cntlegs = cntlegs + 1;
        end
        % Augment speeds counter
        cntspeeds = cntspeeds + 1;
    end
    % Augment trial counter
    cnttrials = cnttrials + 1;
end


end