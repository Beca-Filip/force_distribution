% clear, 
clc, close

P4.ID = 4; % patient ID
P4.Musc_Labels = {'addbrev','addlong','addmagDist','addmagIsch',...
    'addmagMid','addmagProx','glmax1','glmax2','glmax3','glmed1','glmed2',...
    'glmed3','glmin1','glmin2','glmin3','iliacus','psoas','semimem',...
    'semiten','bflh','bfsh','recfem','vasmed','vaslat','vasint',...
    'gaslat','gasmed','tibant','tibpost','perbrev','perlong','pertert',...
    'soleus','edl','fdl'};
P4.nMusc = size(P4.Musc_Labels,2);

% group muscles by their function(s)
[P4.Hip_Flex,P4.Hip_Ext] = muscle_grouping(P4.Musc_Labels,'HipFE');
[P4.Knee_Flex,P4.Knee_Ext] = muscle_grouping(P4.Musc_Labels,'KneeFE');
[P4.Ankle_Df,P4.Ankle_Pf] = muscle_grouping(P4.Musc_Labels,'AnkleFE');

% compute mean and sd of joint stiffness
P4.stiffness.K_L_HipFE = stiffness_processing('HipFE_L');
P4.stiffness.K_L_KneeFE = stiffness_processing('KneeFE_L');
P4.stiffness.K_L_AnkleFE = stiffness_processing('AnkleFE_L');
P4.stiffness.K_R_HipFE = stiffness_processing('HipFE_R');
P4.stiffness.K_R_KneeFE = stiffness_processing('KneeFE_R');
P4.stiffness.K_R_AnkleFE = stiffness_processing('AnkleFE_R');

% compute  contribution to joint stiffness by each muscle
[P4.Contribution_L_HipFE,P4.Stiffness_mus.L_HipFlex,P4.Stiffness_mus.L_HipExt] = stiffness_by_muscle(P4.Hip_Flex,P4.Hip_Ext,'HipFE_L');
[P4.Contribution_L_KneeFE,P4.Stiffness_mus.L_KneeFlex,P4.Stiffness_mus.L_KneeExt] = stiffness_by_muscle(P4.Knee_Flex,P4.Knee_Ext,'KneeFE_L');
[P4.Contribution_L_AnkleFE,P4.Stiffness_mus.L_AnkleDF,P4.Stiffness_mus.L_AnklePF] = stiffness_by_muscle(P4.Ankle_Df,P4.Ankle_Pf,'AnkleFE_L');
[P4.Contribution_R_HipFE,P4.Stiffness_mus.R_HipFlex,P4.Stiffness_mus.R_HipExt] = stiffness_by_muscle(P4.Hip_Flex,P4.Hip_Ext,'HipFE_R');
[P4.Contribution_R_KneeFE,P4.Stiffness_mus.R_KneeFlex,P4.Stiffness_mus.R_KneeExt] = stiffness_by_muscle(P4.Knee_Flex,P4.Knee_Ext,'KneeFE_R');
[P4.Contribution_R_AnkleFE,P4.Stiffness_mus.R_AnkleDF,P4.Stiffness_mus.R_AnklePF] = stiffness_by_muscle(P4.Ankle_Df,P4.Ankle_Pf,'AnkleFE_R');

% process emg
% emg_raw, without calibration for either electromechanical delay or scale
em_delay = 0; scale = 0;
P4.emg_raw = EMG_processing(P4,em_delay,scale);
% emg_scaled, calibrated for scale only, but not electromechanical delay 
em_delay = 0; scale = 1;
P4.emg_scaled = EMG_processing(P4,em_delay,scale);
% emg_delayed, calibrated for electromechanical delay only, but not scale
em_delay = 1; scale = 0;
P4.emg_delayed = EMG_processing(P4,em_delay,scale);
% emg_calibrated, fully calibrated for electromechanical delay and scale
em_delay = 1; scale = 1;
P4.emg_calibrated = EMG_processing(P4,em_delay,scale);

% compute emg-based CCI
[P4.CCI_e_raw,P4.Antagon_Pairs] = CCI_EMG_based(P4.emg_raw,P4);
[P4.CCI_e_scaled,~] = CCI_EMG_based(P4.emg_scaled,P4);
[P4.CCI_e_delayed,~] = CCI_EMG_based(P4.emg_delayed,P4);
[P4.CCI_e_calibrated,~] = CCI_EMG_based(P4.emg_calibrated,P4);

% evalute cor between emg-based CCI and joint stiffness
P4.r_e_raw = evaluate_CCI(P4.CCI_e_raw,P4,'emg');
P4.r_e_scaled = evaluate_CCI(P4.CCI_e_scaled,P4,'emg');
P4.r_e_delayed = evaluate_CCI(P4.CCI_e_delayed,P4,'emg');
P4.r_e_calibrated = evaluate_CCI(P4.CCI_e_calibrated,P4,'emg');

% compute emg-based CCI using CCI formulation 2
P4.CCI_2_e_raw = CCI_EMG_based_2(P4.emg_raw);
P4.CCI_2_e_scaled = CCI_EMG_based_2(P4.emg_scaled);
P4.CCI_2_e_delayed = CCI_EMG_based_2(P4.emg_delayed);
P4.CCI_2_e_calibrated = CCI_EMG_based_2(P4.emg_calibrated);

% evalute cor between emg-based CCI and joint stiffness
P4.r_e_raw_2 = evaluate_CCI(P4.CCI_2_e_raw,P4,'emg');
P4.r_e_scaled_2 = evaluate_CCI(P4.CCI_2_e_scaled,P4,'emg');
P4.r_e_delayed_2 = evaluate_CCI(P4.CCI_2_e_delayed,P4,'emg');
P4.r_e_calibrated_2 = evaluate_CCI(P4.CCI_2_e_calibrated,P4,'emg');

% save data
save('P4.mat','P4')

% functions used by the main script
function [agon,antagon] = muscle_grouping(Musc_Labels,joint)
% this function group muscles by their function, e.g. hip flexor/extensor

% load data from results of EMG-driven model 
load('results_patient4_v22f_optModel5_45678_610sigma_new.mat')

    switch joint
        case 'HipFE'
            MA = HipFEMA;
        case 'KneeFE'
            MA = KneeMA;
        case 'AnkleFE'
            MA = AnkleMA;
    end

% initializing
antagon_idx = [];    
agon_idx = [];

    for i = 1:size(Musc_Labels,2)
        if mean(MA(:,11,i)) < 0 && max(MA(:,11,i)) < 0
            antagon_idx = [antagon_idx, i];
        elseif mean(MA(:,11,i)) > 0 && min(MA(:,11,i)) > 0
            agon_idx = [agon_idx, i];
        end
    end
        
% output value
antagon.idx = antagon_idx;
antagon.labels = Musc_Labels(antagon_idx);
agon.idx = agon_idx;
agon.labels = Musc_Labels(agon_idx);
end

function K = stiffness_processing(joint)
% this function process the joint stiffness data by computing the mean and
% standard deviation of the joint stiffness
trial_idx = 11:20; % trials with self-selected speed of 0.5 m/s

load('Patient4_joint_stiffness_v22f_new.mat')

% extract joint stiffness for the joint of interest
    switch joint
        case 'HipFE_L'
            stiffness = StiffnessLeft.HipFE(:,trial_idx);
        case 'KneeFE_L'
            stiffness = StiffnessLeft.Knee(:,trial_idx);
        case 'AnkleFE_L'
            stiffness = StiffnessLeft.Ankle(:,trial_idx);
        case 'HipFE_R'
            stiffness = StiffnessRight.HipFE(:,trial_idx);
        case 'KneeFE_R'
            stiffness = StiffnessRight.Knee(:,trial_idx);
        case 'AnkleFE_R'
            stiffness = StiffnessRight.Ankle(:,trial_idx);
    end
% store all joint stiffness values
K.value = stiffness;
% compute mean of joint stiffness
K.mean = mean(stiffness,2);
% compute sd of joint stiffness
K.sd = std(stiffness,0,2);
% find the peak joint stiffness value
K.peak = max(max(stiffness));
% find the trough joint stiffness value
K.trough = min(min(stiffness));
end

function [contribution,stiffness_mus_agon,stiffness_mus_antagon] = stiffness_by_muscle(agon,antagon,joint)
% this function computes the mean and sd of contribution to joint stiffness
% by each muscle, as well as the the % of contribution to the overall joint
% stiffness over entire gait cycle by each muscle

trial_idx = 11:20; % trials with self-selected speed of 0.5 m/s
% load individual muscle stiffness data
load('Patient4_joint_stiffness_v22f_individual_mus_new.mat')

% extract stiffness by joint and walking speed
    switch joint
        case 'HipFE_L'
            stiffness_mus_trials = StiffnessLeft_mus.HipFE(:,:,trial_idx);             
        case 'KneeFE_L'
            stiffness_mus_trials = StiffnessLeft_mus.Knee(:,:,trial_idx);
        case 'AnkleFE_L'
            stiffness_mus_trials = StiffnessLeft_mus.Ankle(:,:,trial_idx);
        case 'HipFE_R'
            stiffness_mus_trials = StiffnessRight_mus.HipFE(:,:,trial_idx);
        case 'KneeFE_R'
            stiffness_mus_trials = StiffnessRight_mus.Knee(:,:,trial_idx);
        case 'AnkleFE_R'
            stiffness_mus_trials = StiffnessRight_mus.Ankle(:,:,trial_idx);
    end
    
% extract stiffness by muscle
stiffness_mus_agon = stiffness_mus_trials(agon.idx,:,:);
stiffness_mus_antagon = stiffness_mus_trials(antagon.idx,:,:);
% compute the mean contribution to joint stiffness by each muscle
stiffness_mus_agon_mean = mean(stiffness_mus_agon,3);
stiffness_mus_antagon_mean = mean(stiffness_mus_antagon,3);
% compute the percentage contribution to joint stiffness by each muscle
% area under curve for individual muscles (K vs % gait cycle curve)
stiffness_mus_agon_sum = sum(stiffness_mus_agon_mean,2);
stiffness_mus_antagon_sum = sum(stiffness_mus_antagon_mean,2);
% area under curve for sum of all muscles (K vs % gait cycle curve)
stiffness_mus_sum = sum([stiffness_mus_agon_sum;stiffness_mus_antagon_sum]);
% compute percentage contribution
contribution_agon = stiffness_mus_agon_sum/stiffness_mus_sum;
contribution_antagon = stiffness_mus_antagon_sum/stiffness_mus_sum;
% rank the contributors on both agonist and antagonistic sides
[contribution_agon,idx_agon] = sort(contribution_agon,'descend');
[contribution_antagon,idx_antagon] = sort(contribution_antagon,'descend');
% output results
% fraction of contribution 
contribution.agon_amount = contribution_agon';
contribution.antagon_amount = contribution_antagon';
contribution.agon_idx = agon.idx(idx_agon');
contribution.antagon_idx = antagon.idx(idx_antagon');
% stiffness by individual muscles
stiffness_mus_agon = permute(stiffness_mus_trials(agon.idx(idx_agon'),:,:),[2 3 1]);
stiffness_mus_antagon = permute(stiffness_mus_trials(antagon.idx(idx_antagon'),:,:),[2 3 1]);
end

function emg = EMG_processing(Inpt,em_delay,scale)
% this function process the raw emg and calibrated emg for CCI computation

% indexing
trial_idx_L = 11:20;
trial_idx_R = 61:70;
gait_pts = 21:1:121; % used on raw emg data (141 data points)

% find the various antagonistic pairs for each joint
antagonistic_pairs.L_HipFE = {Inpt.Contribution_L_HipFE.agon_idx(1:end);...
    Inpt.Contribution_L_HipFE.antagon_idx(1:end)};
antagonistic_pairs.L_KneeFE = {Inpt.Contribution_L_KneeFE.agon_idx(1:end);...
    Inpt.Contribution_L_KneeFE.antagon_idx(1:end)};
antagonistic_pairs.L_AnkleFE = {Inpt.Contribution_L_AnkleFE.agon_idx(1:end);...
    Inpt.Contribution_L_AnkleFE.antagon_idx(1:end)};
antagonistic_pairs.R_HipFE = {Inpt.Contribution_R_HipFE.agon_idx(1:end);...
    Inpt.Contribution_R_HipFE.antagon_idx(1:end)};
antagonistic_pairs.R_KneeFE = {Inpt.Contribution_R_KneeFE.agon_idx(1:end);...
    Inpt.Contribution_R_KneeFE.antagon_idx(1:end)};
antagonistic_pairs.R_AnkleFE = {Inpt.Contribution_R_AnkleFE.agon_idx(1:end);...
    Inpt.Contribution_R_AnkleFE.antagon_idx(1:end)};
% this function process emg by the types
% emg_raw - load EMG_raw.mat
% emg_scaled - apply EMGSCaleOpt onto emg_raw
% emg_delayed - descale Excitation by EMGScaleOpt
% emg_calibrated - Excitation
    switch em_delay % either emg_raw or emg_scaled
        case 0
        load('EMG_raw_P4.mat')        
        EMG = permute(EMG,[1 3 2]);     
            switch scale
                case 0 % emg_raw
                    emg.L_HipFlex = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_HipFE{1,:});
                    emg.L_HipExt = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_HipFE{2,:});
                    emg.L_KneeFlex = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_KneeFE{1,:});
                    emg.L_KneeExt = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_KneeFE{2,:});
                    emg.L_AnkleDF = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_AnkleFE{1,:});
                    emg.L_AnklePF = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_AnkleFE{2,:});
                    emg.R_HipFlex = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_HipFE{1,:});
                    emg.R_HipExt = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_HipFE{2,:});
                    emg.R_KneeFlex = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_KneeFE{1,:});
                    emg.R_KneeExt = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_KneeFE{2,:});
                    emg.R_AnkleDF = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_AnkleFE{1,:});
                    emg.R_AnklePF = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_AnkleFE{2,:});
                case 1 % emg_scaled
                    % extract EMG scale, EMGScaleOpt from EMG-driven model
                    load('results_patient4_v22f_optModel5_45678_610sigma_new.mat')
                    % apply scale onto the raw emg signal
                    scale_L_HipFlex = reshape(EMGScaleOpt(antagonistic_pairs.L_HipFE{1,:}),[1 1 size(antagonistic_pairs.L_HipFE{1,:},2)]);
                    emg.L_HipFlex = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_HipFE{1,:}).*scale_L_HipFlex;
                    scale_L_HipExt = reshape(EMGScaleOpt(antagonistic_pairs.L_HipFE{2,:}),[1 1 size(antagonistic_pairs.L_HipFE{2,:},2)]);
                    emg.L_HipExt = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_HipFE{2,:}).*scale_L_HipExt;
                    scale_L_KneeFlex = reshape(EMGScaleOpt(antagonistic_pairs.L_KneeFE{1,:}),[1 1 size(antagonistic_pairs.L_KneeFE{1,:},2)]);
                    emg.L_KneeFlex = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_KneeFE{1,:}).*scale_L_KneeFlex;
                    scale_L_KneeExt = reshape(EMGScaleOpt(antagonistic_pairs.L_KneeFE{2,:}),[1 1 size(antagonistic_pairs.L_KneeFE{2,:},2)]);
                    emg.L_KneeExt = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_KneeFE{2,:}).*scale_L_KneeExt;
                    scale_L_AnkleDF = reshape(EMGScaleOpt(antagonistic_pairs.L_AnkleFE{1,:}),[1 1 size(antagonistic_pairs.L_AnkleFE{1,:},2)]);
                    emg.L_AnkleDF = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_AnkleFE{1,:}).*scale_L_AnkleDF;
                    scale_L_AnklePF = reshape(EMGScaleOpt(antagonistic_pairs.L_AnkleFE{2,:}),[1 1 size(antagonistic_pairs.L_AnkleFE{2,:},2)]);
                    emg.L_AnklePF = EMG(gait_pts,trial_idx_L,antagonistic_pairs.L_AnkleFE{2,:}).*scale_L_AnklePF;                    
                    scale_R_HipFlex = reshape(EMGScaleOpt(antagonistic_pairs.R_HipFE{1,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_HipFE{1,:},2)]);
                    emg.R_HipFlex = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_HipFE{1,:}).*scale_R_HipFlex;                    
                    scale_R_HipExt = reshape(EMGScaleOpt(antagonistic_pairs.R_HipFE{2,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_HipFE{2,:},2)]);
                    emg.R_HipExt = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_HipFE{2,:}).*scale_R_HipExt;                    
                    scale_R_KneeFlex = reshape(EMGScaleOpt(antagonistic_pairs.R_KneeFE{1,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_KneeFE{1,:},2)]);
                    emg.R_KneeFlex = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_KneeFE{1,:}).*scale_R_KneeFlex;                    
                    scale_R_KneeExt = reshape(EMGScaleOpt(antagonistic_pairs.R_KneeFE{2,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_KneeFE{2,:},2)]);
                    emg.R_KneeExt = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_KneeFE{2,:}).*scale_R_KneeExt;                    
                    scale_R_AnkleDF = reshape(EMGScaleOpt(antagonistic_pairs.R_AnkleFE{1,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_AnkleFE{1,:},2)]);
                    emg.R_AnkleDF = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_AnkleFE{1,:}).*scale_R_AnkleDF;                    
                    scale_R_AnklePF = reshape(EMGScaleOpt(antagonistic_pairs.R_AnkleFE{2,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_AnkleFE{2,:},2)]);
                    emg.R_AnklePF = EMG(gait_pts,trial_idx_R,antagonistic_pairs.R_AnkleFE{2,:}).*scale_R_AnklePF;                    
            end
        case 1 % either emg_delayed or emg_calibrated
        load('results_patient4_v22f_optModel5_45678_610sigma_new.mat')
            switch scale
                case 0 % emg_delayed
                    % descale Excitation by EMGScaleOpt
                    scale_L_HipFlex = reshape(EMGScaleOpt(antagonistic_pairs.L_HipFE{1,:}),[1 1 size(antagonistic_pairs.L_HipFE{1,:},2)]);
                    emg.L_HipFlex = Excitations(:,trial_idx_L,antagonistic_pairs.L_HipFE{1,:})./scale_L_HipFlex;                          
                    scale_L_HipExt = reshape(EMGScaleOpt(antagonistic_pairs.L_HipFE{2,:}),[1 1 size(antagonistic_pairs.L_HipFE{2,:},2)]);
                    emg.L_HipExt = Excitations(:,trial_idx_L,antagonistic_pairs.L_HipFE{2,:})./scale_L_HipExt;                    
                    scale_L_KneeFlex = reshape(EMGScaleOpt(antagonistic_pairs.L_KneeFE{1,:}),[1 1 size(antagonistic_pairs.L_KneeFE{1,:},2)]);
                    emg.L_KneeFlex = Excitations(:,trial_idx_L,antagonistic_pairs.L_KneeFE{1,:})./scale_L_KneeFlex;
                    scale_L_KneeExt = reshape(EMGScaleOpt(antagonistic_pairs.L_KneeFE{2,:}),[1 1 size(antagonistic_pairs.L_KneeFE{2,:},2)]);
                    emg.L_KneeExt = Excitations(:,trial_idx_L,antagonistic_pairs.L_KneeFE{2,:})./scale_L_KneeExt;
                    scale_L_AnkleDF = reshape(EMGScaleOpt(antagonistic_pairs.L_AnkleFE{1,:}),[1 1 size(antagonistic_pairs.L_AnkleFE{1,:},2)]);
                    emg.L_AnkleDF = Excitations(:,trial_idx_L,antagonistic_pairs.L_AnkleFE{1,:})./scale_L_AnkleDF;
                    scale_L_AnklePF = reshape(EMGScaleOpt(antagonistic_pairs.L_AnkleFE{2,:}),[1 1 size(antagonistic_pairs.L_AnkleFE{2,:},2)]);
                    emg.L_AnklePF = Excitations(:,trial_idx_L,antagonistic_pairs.L_AnkleFE{2,:})./scale_L_AnklePF;
                    scale_R_HipFlex = reshape(EMGScaleOpt(antagonistic_pairs.R_HipFE{1,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_HipFE{1,:},2)]);
                    emg.R_HipFlex = Excitations(:,trial_idx_R,antagonistic_pairs.R_HipFE{1,:})./scale_R_HipFlex;
                    scale_R_HipExt = reshape(EMGScaleOpt(antagonistic_pairs.R_HipFE{2,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_HipFE{2,:},2)]);
                    emg.R_HipExt = Excitations(:,trial_idx_R,antagonistic_pairs.R_HipFE{2,:})./scale_R_HipExt;
                    scale_R_KneeFlex = reshape(EMGScaleOpt(antagonistic_pairs.R_KneeFE{1,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_KneeFE{1,:},2)]);
                    emg.R_KneeFlex = Excitations(:,trial_idx_R,antagonistic_pairs.R_KneeFE{1,:})./scale_R_KneeFlex;
                    scale_R_KneeExt = reshape(EMGScaleOpt(antagonistic_pairs.R_KneeFE{2,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_KneeFE{2,:},2)]);
                    emg.R_KneeExt = Excitations(:,trial_idx_R,antagonistic_pairs.R_KneeFE{2,:})./scale_R_KneeExt;
                    scale_R_AnkleDF = reshape(EMGScaleOpt(antagonistic_pairs.R_AnkleFE{1,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_AnkleFE{1,:},2)]);
                    emg.R_AnkleDF = Excitations(:,trial_idx_R,antagonistic_pairs.R_AnkleFE{1,:})./scale_R_AnkleDF;
                    scale_R_AnklePF = reshape(EMGScaleOpt(antagonistic_pairs.R_AnkleFE{2,:}+Inpt.nMusc),[1 1 size(antagonistic_pairs.R_AnkleFE{2,:},2)]);
                    emg.R_AnklePF = Excitations(:,trial_idx_R,antagonistic_pairs.R_AnkleFE{2,:})./scale_R_AnklePF;
                case 1 % emg_calibrated
                    emg.L_HipFlex = Excitations(:,trial_idx_L,antagonistic_pairs.L_HipFE{1,:});
                    emg.L_HipExt = Excitations(:,trial_idx_L,antagonistic_pairs.L_HipFE{2,:});
                    emg.L_KneeFlex = Excitations(:,trial_idx_L,antagonistic_pairs.L_KneeFE{1,:});
                    emg.L_KneeExt = Excitations(:,trial_idx_L,antagonistic_pairs.L_KneeFE{2,:});
                    emg.L_AnkleDF = Excitations(:,trial_idx_L,antagonistic_pairs.L_AnkleFE{1,:});
                    emg.L_AnklePF = Excitations(:,trial_idx_L,antagonistic_pairs.L_AnkleFE{2,:});
                    emg.R_HipFlex = Excitations(:,trial_idx_R,antagonistic_pairs.R_HipFE{1,:});
                    emg.R_HipExt = Excitations(:,trial_idx_R,antagonistic_pairs.R_HipFE{2,:});
                    emg.R_KneeFlex = Excitations(:,trial_idx_R,antagonistic_pairs.R_KneeFE{1,:});
                    emg.R_KneeExt = Excitations(:,trial_idx_R,antagonistic_pairs.R_KneeFE{2,:});
                    emg.R_AnkleDF = Excitations(:,trial_idx_R,antagonistic_pairs.R_AnkleFE{1,:});
                    emg.R_AnklePF = Excitations(:,trial_idx_R,antagonistic_pairs.R_AnkleFE{2,:});                   
            end
    end
end

function [CCI_e,Antagon_Pairs] = CCI_EMG_based(emg,Inpt)
% this function computes CCI based on emg data of antagonistic pairs
% number of muscles in each antagonistic group used for CCI computation
n_HipFlex = size(emg.L_HipFlex,3);
n_HipExt = size(emg.L_HipExt,3);
n_KneeFlex = size(emg.L_KneeFlex,3);
n_KneeExt = size(emg.L_KneeExt,3);
n_AnkleDF = size(emg.L_AnkleDF,3);
n_AnklePF = size(emg.L_AnklePF,3);
% form combinations of different antagonistic pairs
combo_HipFE = combvec([1:n_HipFlex],[1:n_HipExt]);
combo_KneeFE = combvec([1:n_KneeFlex],[1:n_KneeExt]);
combo_AnkleFE = combvec([1:n_AnkleDF],[1:n_AnklePF]);
% get the global muscle index and label for each combination index
Antagon_Pairs.HipFE_idx = [Inpt.Contribution_L_HipFE.agon_idx(combo_HipFE(1,:));...
    Inpt.Contribution_L_HipFE.antagon_idx(combo_HipFE(2,:))];
Antagon_Pairs.KneeFE_idx = [Inpt.Contribution_L_KneeFE.agon_idx(combo_KneeFE(1,:));...
    Inpt.Contribution_L_KneeFE.antagon_idx(combo_KneeFE(2,:))];
Antagon_Pairs.AnkleFE_idx = [Inpt.Contribution_L_AnkleFE.agon_idx(combo_AnkleFE(1,:));...
    Inpt.Contribution_L_AnkleFE.antagon_idx(combo_AnkleFE(2,:))];
% label
Antagon_Pairs.HipFE_label = Inpt.Musc_Labels(Antagon_Pairs.HipFE_idx);
Antagon_Pairs.KneeFE_label = Inpt.Musc_Labels(Antagon_Pairs.KneeFE_idx);
Antagon_Pairs.AnkleFE_label = Inpt.Musc_Labels(Antagon_Pairs.AnkleFE_idx);

    for i = 1:size(combo_HipFE,2) % each agon-antagon combination, HipFE
        for j = 1:size(emg.L_HipFlex,2) % number of trials
            CCI_e.L_HipFE(:,j,i) = CCI_Formulation1(emg.L_HipFlex(:,j,combo_HipFE(1,i)),emg.L_HipExt(:,j,combo_HipFE(2,i)));
            CCI_e.R_HipFE(:,j,i) = CCI_Formulation1(emg.R_HipFlex(:,j,combo_HipFE(1,i)),emg.R_HipExt(:,j,combo_HipFE(2,i)));            
        end
    end
    for i = 1:size(combo_KneeFE,2) % each agon-antagon combination, KneeFE
        for j = 1:size(emg.L_KneeFlex,2) % number of trials
            CCI_e.L_KneeFE(:,j,i) = CCI_Formulation1(emg.L_KneeFlex(:,j,combo_KneeFE(1,i)),emg.L_KneeExt(:,j,combo_KneeFE(2,i)));
            CCI_e.R_KneeFE(:,j,i) = CCI_Formulation1(emg.R_KneeFlex(:,j,combo_KneeFE(1,i)),emg.R_KneeExt(:,j,combo_KneeFE(2,i)));
        end
    end
    for i = 1:size(combo_AnkleFE,2) % each possible agon-antagon combination
        for j = 1:size(emg.L_AnkleDF,2) % number of trials
            CCI_e.L_AnkleFE(:,j,i) = CCI_Formulation1(emg.L_AnkleDF(:,j,combo_AnkleFE(1,i)),emg.L_AnklePF(:,j,combo_AnkleFE(2,i)));
            CCI_e.R_AnkleFE(:,j,i) = CCI_Formulation1(emg.R_AnkleDF(:,j,combo_AnkleFE(1,i)),emg.R_AnklePF(:,j,combo_AnkleFE(2,i)));
        end
    end    
end

function CCI_e = CCI_EMG_based_2(emg)

n_HipFlex = size(emg.L_HipFlex,3);
n_HipExt = size(emg.L_HipExt,3);
n_KneeFlex = size(emg.L_KneeFlex,3);
n_KneeExt = size(emg.L_KneeExt,3);
n_AnkleDF = size(emg.L_AnkleDF,3);
n_AnklePF = size(emg.L_AnklePF,3);

% form combinations of different antagonistic pairs
combo_HipFE = combvec([1:n_HipFlex],[1:n_HipExt]);
combo_KneeFE = combvec([1:n_KneeFlex],[1:n_KneeExt]);
combo_AnkleFE = combvec([1:n_AnkleDF],[1:n_AnklePF]);

    for i = 1:size(combo_HipFE,2) % each agon-antagon combination, HipFE
        for j = 1:size(emg.L_HipFlex,2)
            CCI_e.L_HipFE(:,j,i) = CCI_Formulation2(emg.L_HipFlex(:,j,combo_HipFE(1,i)),emg.L_HipExt(:,j,combo_HipFE(2,i)));
            CCI_e.R_HipFE(:,j,i) = CCI_Formulation2(emg.R_HipFlex(:,j,combo_HipFE(1,i)),emg.R_HipExt(:,j,combo_HipFE(2,i)));                    
        end
    end
    for i = 1:size(combo_KneeFE,2) % each agon-antagon combination, KneeFE
        for j = 1:size(emg.L_KneeFlex,2) % number of trials
            CCI_e.L_KneeFE(:,j,i) = CCI_Formulation2(emg.L_KneeFlex(:,j,combo_KneeFE(1,i)),emg.L_KneeExt(:,j,combo_KneeFE(2,i)));
            CCI_e.R_KneeFE(:,j,i) = CCI_Formulation2(emg.R_KneeFlex(:,j,combo_KneeFE(1,i)),emg.R_KneeExt(:,j,combo_KneeFE(2,i)));
        end
    end
    for i = 1:size(combo_AnkleFE,2) % each possible agon-antagon combination           
        for j = 1:size(emg.L_AnkleDF,2) % number of trials
            CCI_e.L_AnkleFE(:,j,i) = CCI_Formulation2(emg.L_AnkleDF(:,j,combo_AnkleFE(1,i)),emg.L_AnklePF(:,j,combo_AnkleFE(2,i)));
            CCI_e.R_AnkleFE(:,j,i) = CCI_Formulation2(emg.R_AnkleDF(:,j,combo_AnkleFE(1,i)),emg.R_AnklePF(:,j,combo_AnkleFE(2,i)));
        end
    end
end

function CCI = CCI_Formulation1(Inpt1,Inpt2)
% this function uses the Rudolph formulation to compute CCI
% initializing
CCI = zeros(size(Inpt1));
H_Inpt = zeros(size(CCI)); % input with higher value
L_Inpt = zeros(size(CCI)); % input with lower value
% obtain the absolute value of each input
Inpt1 = abs(Inpt1);
Inpt2 = abs(Inpt2);
% finding where Inpt1 is greater than Inpt2 and vice versa
idx1 = find(Inpt1>Inpt2);
idx2 = find(Inpt1<Inpt2);
idx3 = find(Inpt1==Inpt2);
H_Inpt(idx1) = Inpt1(idx1);
L_Inpt(idx2) = Inpt1(idx2);
H_Inpt(idx2) = Inpt2(idx2);
L_Inpt(idx1) = Inpt2(idx1);
H_Inpt(idx3) = Inpt1(idx3);
L_Inpt(idx3) = Inpt2(idx3);
% compute CCI using Rudolph formulation
CCI = L_Inpt/H_Inpt*(L_Inpt+H_Inpt);
end

function CCI = CCI_Formulation2(Inpt1,Inpt2)
% this function uses the Falconer and Winter formulation to compute CCI
% initializing
CCI = zeros(size(Inpt1));
H_Inpt = zeros(size(CCI)); % input with higher value
L_Inpt = zeros(size(CCI)); % input with lower value
% obtain the absolute value of each input
Inpt1 = abs(Inpt1);
Inpt2 = abs(Inpt2);
% finding where Inpt1 is greater than Inpt2 and vice versa
idx1 = find(Inpt1>Inpt2);
idx2 = find(Inpt1<Inpt2);
idx3 = find(Inpt1==Inpt2);
H_Inpt(idx1) = Inpt1(idx1);
L_Inpt(idx2) = Inpt1(idx2);
H_Inpt(idx2) = Inpt2(idx2);
L_Inpt(idx1) = Inpt2(idx1);
H_Inpt(idx3) = Inpt1(idx3);
L_Inpt(idx3) = Inpt2(idx3);
% compute CCI using Rudolph formulation
CCI = 2*L_Inpt./(L_Inpt+H_Inpt);
end

function r = evaluate_CCI(CCI,Inpt,CCI_type)
% this function computes the correlation between CCI and joint stiffness
    switch CCI_type
        case 'emg'
            for i = 1:size(CCI.L_HipFE,3) % antagonistic pairs, Hip FE
                for j = 1:size(CCI.L_HipFE,2) % number of trials
                    R_L_HipFE = corrcoef(CCI.L_HipFE(:,j,i),Inpt.stiffness.K_L_HipFE.value(:,j));
                    r.L_HipFE(j,i) = R_L_HipFE(1,2);
                    R_R_HipFE = corrcoef(CCI.R_HipFE(:,j,i),Inpt.stiffness.K_R_HipFE.value(:,j));
                    r.R_HipFE(j,i) = R_R_HipFE(1,2);                        
                end
            end
            for i = 1:size(CCI.L_KneeFE,3) % antagonistic pairs, KneeFE
                for j = 1:size(CCI.L_KneeFE,2) % number of trials
                    R_L_KneeFE = corrcoef(CCI.L_KneeFE(:,j,i),Inpt.stiffness.K_L_KneeFE.value(:,j));
                    r.L_KneeFE(j,i) = R_L_KneeFE(1,2);
                    R_R_KneeFE = corrcoef(CCI.R_KneeFE(:,j,i),Inpt.stiffness.K_R_KneeFE.value(:,j));
                    r.R_KneeFE(j,i) = R_R_KneeFE(1,2);
                end
            end
            for i = 1:size(CCI.L_AnkleFE,3) % antagonistic pairs, AnkleFE
                for j = 1:size(CCI.L_AnkleFE,2) % number of trials
                    R_L_AnkleFE = corrcoef(CCI.L_AnkleFE(:,j,i),Inpt.stiffness.K_L_AnkleFE.value(:,j));
                    r.L_AnkleFE(j,i) = R_L_AnkleFE(1,2);
                    R_R_AnkleFE = corrcoef(CCI.R_AnkleFE(:,j,i),Inpt.stiffness.K_R_AnkleFE.value(:,j));
                    r.R_AnkleFE(j,i) = R_R_AnkleFE(1,2);
                end
            end
        case 'moment'
            for j = 1:size(CCI.L_HipFE,2) % number of trials
                R_L_HipFE = corrcoef(CCI.L_HipFE(:,j),Inpt.stiffness.K_L_HipFE.value(:,j));
                r.L_HipFE(j,1) = R_L_HipFE(1,2);
                R_L_KneeFE = corrcoef(CCI.L_KneeFE(:,j),Inpt.stiffness.K_L_KneeFE.value(:,j));
                r.L_KneeFE(j,1) = R_L_KneeFE(1,2);
                R_L_AnkleFE = corrcoef(CCI.L_AnkleFE(:,j),Inpt.stiffness.K_L_AnkleFE.value(:,j));
                r.L_AnkleFE(j,1) = R_L_AnkleFE(1,2);
                R_R_HipFE = corrcoef(CCI.R_HipFE(:,j),Inpt.stiffness.K_R_HipFE.value(:,j));
                r.R_HipFE(j,1) = R_R_HipFE(1,2);
                R_R_KneeFE = corrcoef(CCI.R_KneeFE(:,j),Inpt.stiffness.K_R_KneeFE.value(:,j));
                r.R_KneeFE(j,1) = R_R_KneeFE(1,2);
                R_R_AnkleFE = corrcoef(CCI.R_AnkleFE(:,j),Inpt.stiffness.K_R_AnkleFE.value(:,j));
                r.R_AnkleFE(j,1) = R_R_AnkleFE(1,2);
            end
    end
end

