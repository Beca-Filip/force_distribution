%%  ***Stiffness calculation function***
%
%  This function computes joint muscular stiffness absed on current
%  muscle fiber length l and velocity v and activation level of each muscle.
%  Compute_Joint_Stiffness_v3:
%  add calculation of individual muscle contribution to the joint stiffness
%  choose consistent method and parameters for muscle stiffness calculation


function  [JointStiffness,JointStiffnessMus,dfdl] = Compute_Joint_Stiffness_v4(a,l,v,r,drdq,F,MusParams)

% extract variables from MusParams
nMusc = MusParams.nMusc; % number of muscles
nframesAll = MusParams.nframesAll; % number of all time frames
lmoOpt = MusParams.lmoOpt ; % lmo
ltsOpt = MusParams.ltsOpt ; % lmo
Fmax = MusParams.Fmax; % Fmax
alpha = MusParams.alpha; % alpha
VmaxFactor = MusParams.VmaxFactor; %VmaxFactor = 10

% preallocate memory for joint stiffness results
JointStiffness = zeros(1,nframesAll);
JointStiffnessMus = zeros(nMusc, nframesAll);
dfdl = zeros(nMusc, nframesAll);


for j = 1:nMusc
%     for i = 1:nframesAll
        
        % Construct 'Params' structure for dfdl calculation
        % Extract instaneous values for computation of joint stiffness
        
        Params.Fmax_sym = Fmax(j);
        Params.lmo_sym = lmoOpt(j);
        Params.lts_sym = ltsOpt(j);
        Params.alphaAngle = alpha(j);
        Params.VmaxFactor_sym = VmaxFactor;
        Params.Lmt_sym = l(:,j);
        Params.Vmt_sym = v(:,j);
        Params.a_sym = a(:,j);
        
        % Calculate dfdl values
        dfdl(j,:) = dfdlcal_v2(Params);
        
        %moment arms
        r_ins = r(:,j); r_ins = r_ins';
        %muscle force
        Force_ins = F(:,j); Force_ins = Force_ins';
        %drdq
        drdq_ins = drdq(:,j); drdq_ins = drdq_ins';
%         for i = 1:nframesAll        
            % Calculate joint stiffness
            JointStiffnessMus(j,:) = -(drdq_ins .* Force_ins - r_ins.^2 .* dfdl(j,:));
            JointStiffness = JointStiffness + JointStiffnessMus(j,:);
%         end
%     end
end
end

