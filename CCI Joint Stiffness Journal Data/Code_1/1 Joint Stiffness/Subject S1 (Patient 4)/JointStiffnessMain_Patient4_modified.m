% Calculation of joint stiffness
% clearvars -except leg_list speed_list leg_list trial speed cycle; close all;

%% Load optimization results
LoadFileName = 'results_patient4_v22f_optModel5_45678_610sigma_new.mat'; % all results for the simulation get saved here
load(LoadFileName);

%% Initial setting
SampleStep = 1; % resamples data to ever SampleStep points using 1:SampleStep:end. Must be 1, 2, 5 or 10.
nptsLong = (141-1)/SampleStep+1; % number of sample points in gait cycle with padded beginnings and ends
nptsShort = (101-1)/SampleStep+1; % number of sample points in gait cycle without padding
numTrialsPerSpeed = 10; % number of trials used for each speed
Speeds = [0.4 0.5 0.6 0.7 0.8]; % which speeds will be loaded
Leg = 'both'; %  Determine number of gait speeds to be used in later optimizations
nMusc = 35; % Specify number of muscles and EMG signals for each muscle
SaveStiffness = 0; % 1; % save stiffness results
SaveFileName1 = 'Patient4_joint_stiffness_v22f_new.mat';
SaveFileName2 = 'Patient4_joint_stiffness_v22f_individual_mus_new.mat';

%% Initialize Data
% Determine number of gait speeds to be used in later optimizations
if strcmpi(Leg,'both')
    nSpeeds = 2*length(Speeds);
elseif strcmpi(Leg,'left')||strcmpi(Leg,'right')
    nSpeeds = length(Speeds);
end
nTrials = nSpeeds*numTrialsPerSpeed;

% Reform the data
Lmt = reshape(Lmt, nptsShort*nTrials, nMusc);
Vmt= reshape(permute(Vmt, [1 3 2]), nptsShort*nTrials, nMusc);
HipAAMA = reshape(HipAAMA, nptsShort*nTrials, nMusc);
HipFEMA = reshape(HipFEMA, nptsShort*nTrials, nMusc);
KneeMA = reshape(KneeMA, nptsShort*nTrials, nMusc);
AnkleMA = reshape(AnkleMA, nptsShort*nTrials, nMusc);
SubtalarMA =  reshape(SubtalarMA, nptsShort*nTrials, nMusc);
MuscleForce = reshape(muscleForce, nptsShort*nTrials, nMusc);
JAngles = reshape(permute(JAngles, [1 3 2]), nptsShort*nTrials, 6);
JVels = reshape(permute(JVels, [1 3 2]), nptsShort*nTrials, 6);

%% Assign muscles to the actuators of joint movement
% define muscle reference matrix that defines which joint(s) a muscle
% actuates
MuscRef = zeros(1,nMusc);
MuscRef(1:17) = 1; % Hip FE and AA
MuscRef([18:20 22]) = 2; % Hip FE and AA and Knee FE
MuscRef([21 23:25]) = 3; % Knee FE
MuscRef([26:27]) = 4; % Knee FE, Ankle and Subtalar
MuscRef([28:35]) = 5; % Ankle and Subtalar

%% Calculation of drdq
MatDrDq = buildDrDqMatrices(JAngles);
MuscRef = MuscRef;
k = 1;
nframesAll = nptsShort*nTrials;
Ref1 = 1:nframesAll;
Ref2 = nframesAll+1:2*nframesAll;
Ref3 = 2*nframesAll+1:3*nframesAll;

%Preallocate memory
HipFE_drdq = zeros(nframesAll, nMusc);
HipAA_drdq = zeros(nframesAll, nMusc);
Knee_drdq = zeros(nframesAll, nMusc);
Ankle_drdq = zeros(nframesAll, nMusc);
Subtalar_drdq = zeros(nframesAll, nMusc);

% Calculation of drdq; only 5 type of muscles
for i = 1:nMusc
    if MuscRef(i) == 1
        Vec = MatDrDq {1}*coefsNew(k:k+18);
        HipFE_drdq(:,i) = Vec(Ref1);
        HipAA_drdq(:,i) = Vec(Ref2);
        k = k+19;
    elseif MuscRef(i) == 2
        Vec = MatDrDq{2}*coefsNew(k:k+21);
        HipFE_drdq(:,i) = Vec(Ref1);
        HipAA_drdq(:,i) = Vec(Ref2);
        Knee_drdq(:,i) = Vec(Ref3);
        k = k+22;
    elseif MuscRef(i) == 3
        Vec = MatDrDq{3}*coefsNew(k:k+3);
        Knee_drdq(:,i) = Vec(Ref1);
        k = k+4;
    elseif MuscRef(i) == 4
        Vec = MatDrDq{4}*coefsNew(k:k+12);
        Knee_drdq(:,i) = Vec(Ref1);
        Ankle_drdq(:,i) = Vec(Ref2);
        Subtalar_drdq(:,i) = Vec(Ref3);
        k = k+13;
    elseif MuscRef(i) == 5
        Vec = MatDrDq{5}*coefsNew(k:k+9);
        Ankle_drdq(:,i) = Vec(Ref1);
        Subtalar_drdq(:,i) = Vec(Ref2);
        k = k+10;
    end
end

%% Construct 'muscleParams'
MusParams.nMusc = nMusc; % number of muscles
MusParams.nframesAll = nframesAll; % number of all time frames
MusParams.lmoOpt = lmoOpt; % lmo
MusParams.ltsOpt = ltsOpt; % lmo
MusParams.Fmax = Fmax; % Fmax
MusParams.alpha = alpha; % alpha
MusParams.VmaxFactor = optParams.VmaxFactor; % VmaxFactor

%% Calculate joint stiffness
[Jstiffness_HipFE,    Jstiffness_HipFE_mus,    HipFE_dFdl] = Compute_Joint_Stiffness_v4(a,Lmt,Vmt,HipFEMA,HipFE_drdq,MuscleForce,MusParams);
[Jstiffness_HipAA,    Jstiffness_HipAA_mus,    HipAA_dFdl] = Compute_Joint_Stiffness_v4(a,Lmt,Vmt,HipAAMA,HipAA_drdq,MuscleForce,MusParams);
[Jstiffness_Knee,     Jstiffness_Knee_mus,     Knee_dFdl] = Compute_Joint_Stiffness_v4(a,Lmt,Vmt,KneeMA,Knee_drdq,MuscleForce,MusParams);
[Jstiffness_Ankle,    Jstiffness_Ankle_mus,    Ankle_dFdl] = Compute_Joint_Stiffness_v4(a,Lmt,Vmt,AnkleMA,Ankle_drdq,MuscleForce,MusParams);
[Jstiffness_Subtalar, Jstiffness_Subtalar_mus,Subtalar_dFdl] = Compute_Joint_Stiffness_v4(a,Lmt,Vmt,SubtalarMA,Subtalar_drdq,MuscleForce,MusParams);


% Organize data by trial and side
% Left
for i = 1 : nTrials/2
    StiffnessLeft.HipAA(:,i)= Jstiffness_HipAA((i-1)*nptsShort+1: i*nptsShort);
    StiffnessLeft.HipFE(:,i) = Jstiffness_HipFE((i-1)*nptsShort+1: i*nptsShort);
    StiffnessLeft.Knee(:,i) = Jstiffness_Knee((i-1)*nptsShort+1: i*nptsShort);
    StiffnessLeft.Ankle(:,i) = Jstiffness_Ankle((i-1)*nptsShort+1: i*nptsShort);
    StiffnessLeft.Subtalar(:,i) = Jstiffness_Subtalar((i-1)*nptsShort+1: i*nptsShort);
end
for i = 1 : nTrials/2
    StiffnessLeft_mus.HipAA(:,:,i)= Jstiffness_HipAA_mus(:,(i-1)*nptsShort+1: i*nptsShort);
    StiffnessLeft_mus.HipFE(:,:,i) = Jstiffness_HipFE_mus(:,(i-1)*nptsShort+1: i*nptsShort);
    StiffnessLeft_mus.Knee(:,:,i) = Jstiffness_Knee_mus(:,(i-1)*nptsShort+1: i*nptsShort);
    StiffnessLeft_mus.Ankle(:,:,i) = Jstiffness_Ankle_mus(:,(i-1)*nptsShort+1: i*nptsShort);
    StiffnessLeft_mus.Subtalar(:,:,i) = Jstiffness_Subtalar_mus(:,(i-1)*nptsShort+1: i*nptsShort);
end
% Right
for i =  nTrials/2+1 : nTrials
    StiffnessRight.HipAA(:,i-nTrials/2)= Jstiffness_HipAA((i-1)*nptsShort+1: i*nptsShort);
    StiffnessRight.HipFE(:,i-nTrials/2) = Jstiffness_HipFE((i-1)*nptsShort+1: i*nptsShort);
    StiffnessRight.Knee(:,i-nTrials/2) = Jstiffness_Knee((i-1)*nptsShort+1: i*nptsShort);
    StiffnessRight.Ankle(:,i-nTrials/2) = Jstiffness_Ankle((i-1)*nptsShort+1: i*nptsShort);
    StiffnessRight.Subtalar(:,i-nTrials/2) = Jstiffness_Subtalar((i-1)*nptsShort+1: i*nptsShort);
end
for i =  nTrials/2+1 : nTrials
    StiffnessRight_mus.HipAA(:,:,i-nTrials/2)= Jstiffness_HipAA_mus(:,(i-1)*nptsShort+1: i*nptsShort);
    StiffnessRight_mus.HipFE(:,:,i-nTrials/2) = Jstiffness_HipFE_mus(:,(i-1)*nptsShort+1: i*nptsShort);
    StiffnessRight_mus.Knee(:,:,i-nTrials/2) = Jstiffness_Knee_mus(:,(i-1)*nptsShort+1: i*nptsShort);
    StiffnessRight_mus.Ankle(:,:,i-nTrials/2) = Jstiffness_Ankle_mus(:,(i-1)*nptsShort+1: i*nptsShort);
    StiffnessRight_mus.Subtalar(:,:,i-nTrials/2) = Jstiffness_Subtalar_mus(:,(i-1)*nptsShort+1: i*nptsShort);
end

% save data
if SaveStiffness ==1
    save(SaveFileName1,'StiffnessLeft', 'StiffnessRight');
    save(SaveFileName2,'StiffnessLeft_mus','StiffnessRight_mus');
end

 

%% functions
function [Mat_drdq] = buildDrDqMatrices(JAngles)
nFrames = size(JAngles,1);
onesCol = ones(nFrames,1);
zerosMat = zeros(nFrames,3);

% matrices for determining drdq
HipFE_drdq = -[2*onesCol 6*JAngles(:,1)];
HipAA_drdq = -[2*onesCol 6*JAngles(:,2)];
InteractionFEAA_drdq_FE = -[2*JAngles(:,2) 0*onesCol];
InteractionFEAA_drdq_AA = -[0*onesCol 2*JAngles(:,1)];
InteractionFEIE_drdq_FE = -[2*JAngles(:,6) 0*onesCol];
InteractionIEAA_drdq_AA = -[ 0*onesCol 2*JAngles(:,6)];
Knee_drdq = -[2*onesCol 6*JAngles(:,3)];
Ankle_drdq = -[2*onesCol 6*JAngles(:,4)];
Subtalar_drdq = -[2*onesCol 6*JAngles(:,5)];
InteractionAnkleSub_drdq_ankle = -[2*JAngles(:,5) 0*onesCol];
InteractionAnkleSub_drdq_Sub = -[0*onesCol 2*JAngles(:,4)];

%HipFE/AA(2 row, and in each row, the direvative w.r.t FE or AA)
Mat_drdq{1} = [0*onesCol 0*onesCol HipFE_drdq zerosMat zerosMat 0*onesCol InteractionFEAA_drdq_FE zerosMat 0*onesCol InteractionFEIE_drdq_FE;...
    0*onesCol zerosMat 0*onesCol HipAA_drdq zerosMat 0*onesCol InteractionFEAA_drdq_AA 0*onesCol InteractionIEAA_drdq_AA zerosMat];
% %HipFE/AA/Knee
Mat_drdq{2} = [0*onesCol 0*onesCol HipFE_drdq zerosMat zerosMat zerosMat 0*onesCol InteractionFEAA_drdq_FE zerosMat 0*onesCol InteractionFEIE_drdq_FE;...
    0*onesCol zerosMat 0*onesCol HipAA_drdq zerosMat zerosMat 0*onesCol InteractionFEAA_drdq_AA 0*onesCol InteractionIEAA_drdq_AA zerosMat;...
    0*onesCol zerosMat zerosMat zerosMat 0*onesCol Knee_drdq zerosMat zerosMat zerosMat];
%Knee
Mat_drdq{3} = [0*onesCol 0*onesCol Knee_drdq];
%Knee/Ankle/Sub
Mat_drdq{4} =  [0*onesCol 0*onesCol Knee_drdq zerosMat zerosMat zerosMat;...
    0*onesCol zerosMat 0*onesCol Ankle_drdq zerosMat 0*onesCol InteractionAnkleSub_drdq_ankle;...
    0*onesCol zerosMat zerosMat 0*onesCol Subtalar_drdq 0*onesCol InteractionAnkleSub_drdq_Sub];
% %Ankle/Sub
Mat_drdq{5} =  [0*onesCol 0*onesCol Ankle_drdq zerosMat 0*onesCol InteractionAnkleSub_drdq_ankle;...
    0*onesCol zerosMat 0*onesCol Subtalar_drdq 0*onesCol InteractionAnkleSub_drdq_Sub];
end





