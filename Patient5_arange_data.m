close all;
clear all;
clc

%% Cycle

% In original script
% trial = 5; % 1 to 10 (per speed)
% leg = 1; % 1 contralateral (left-side) leg / 2 hemiparetic (right-side) leg
% speed = 5; % 2 to 6 (0.25 to 0.65 m/s).
% cycle = (leg-1)*50 + (speed-2)*10 + trial; % Corresponding cycle in results

% Data number of instances
ntrials = 10;
nlegs = 2;
nspeeds = 5;
trial_list = 1:10;
leg_list = 1:2;
speed_list = 2:6;

% Dimensionalities
n = 34;     % Dimensionality of forces
ne = 5;     % Dimensionality of moments (num. equality constraints)
T = 101;    % Number of time samples

% Create a struct to store data
data.description = strvcat(" ",...
                          "All fields have their first two slices equal to data matrix form.",...
                          "-Fields with 5 slices total are constant w.r.t. time and slices 3,4,5 correspond to trial,speed,leg respectively.", ...
                          "-Fields with 6 slices total are non-constant in time and slice 3 correspond to time-samples (of which there are 101) and slices 4,5,6 correspond to trial,speed,leg respectively.", ...
                          " ", ...
                          "Fields J_max and J_min containt cost function normalization data and have size 1x17.");

% Create tensor fields to store data
% Time-variable data
% Dimensionality same as force
data.f = zeros(n, 1, T, ntrials, nspeeds, nlegs);
data.fmin = zeros(n, 1, T, ntrials, nspeeds, nlegs);
data.fmax = zeros(n, 1, T, ntrials, nspeeds, nlegs);
data.lmt = zeros(n, 1, T, ntrials, nspeeds, nlegs);
data.vmt = zeros(n, 1, T, ntrials, nspeeds, nlegs);
data.M = zeros(n, 1, T, ntrials, nspeeds, nlegs);
data.fpassive = zeros(n, 1, T, ntrials, nspeeds, nlegs);
% Dimensionality same as moment
data.A_orig = zeros(ne, n, T, ntrials, nspeeds, nlegs);
data.A = zeros(ne, n, T, ntrials, nspeeds, nlegs);
data.b = zeros(ne, 1, T, ntrials, nspeeds, nlegs);
data.drdq = zeros(ne, n, T, ntrials, nspeeds, nlegs);
data.Kref = zeros(ne, 1, T, ntrials, nspeeds, nlegs);
data.K = zeros(ne, 1, T, ntrials, nspeeds, nlegs);

% Time-constant data
% Dimensionality same as force
data.pcsa = zeros(n, 1, ntrials, nspeeds, nlegs);
data.f0 = zeros(n, 1, ntrials, nspeeds, nlegs);
data.mass = zeros(n, 1, ntrials, nspeeds, nlegs);
data.r = zeros(n, 1, ntrials, nspeeds, nlegs);
data.lmo = zeros(n, 1, ntrials, nspeeds, nlegs);
data.lts = zeros(n, 1, ntrials, nspeeds, nlegs);
data.penation = zeros(n, 1, ntrials, nspeeds, nlegs);

% Dimensionality one
data.foot_off_index_m = zeros(1, 1, ntrials, nspeeds, nlegs); % Computed by moment
data.foot_off_index_a = zeros(1, 1, ntrials, nspeeds, nlegs); % Computed by angle

% Cost function normalization
data.J_max = -inf*ones(1, 17);
data.J_min = inf*ones(1, 17);
%%
for leg = leg_list
    for speed = speed_list
        for trial = trial_list
            %% Calculate corresponding cycle
            cycle = (leg-1)*50 + (speed-2)*10 + trial; % Corresponding cycle in results
            %% Paths to data and code
            addpath('CCI Joint Stiffness Journal Data\Code_1\1 Joint Stiffness\Subject S2 (Patient 5)')
            addpath('CCI Joint Stiffness Journal Data\Processed Data\Patient 5')
            %% Run code to compute joint stiffness
            JointStiffnessMain_Patient5_modified
            % All derivatives are computed analytically (lever arms and muscle lenghts 
            % are approximated using cubic polynomial equations)
            % This code first loads results_patient4_v22f_optModel5_45678_610sigma_new.mat
            % This code reshapes HipAAMA, HipFEMA, KneeMA, AnkleMA, SubtalarMA, MuscleForce
            % All variables are cleared except the ones needed to compute joint stiffness
            clearvars -except cycle ...
                HipFE_drdq HipAA_drdq Knee_drdq Ankle_drdq Subtalar_drdq ...
                HipFE_dFdl JAngles alpha...
                n ne T...
                Jstiffness_HipFE Jstiffness_HipAA Jstiffness_Knee ...
                Jstiffness_Ankle Jstiffness_Subtalar ...
                leg_list speed_list trial_list leg speed trial data ...                
            % Same results in HipFE_dFdl, HipAA_dFdl, Knee_dFdl, Ankle_dFdl, or Subtalar_dFdl
            %% Load gait analysis data
            vars = who('-file','results_patient5_v22f_optModel5_23456_610sigma_new.mat');
            vars_to_load = setdiff(vars, 'JAngles');
            load('results_patient5_v22f_optModel5_23456_610sigma_new.mat', vars_to_load{:});

            % Joint angles
            theta = JAngles((cycle-1)*101+1:cycle*101,[1:5]);
            
            % Joint moments
            b(:,1) = IDloads(:,cycle,1); % Hip Flexion / Extension
            b(:,2) = IDloads(:,cycle,2); % Hip Adduction / Abduction
            b(:,3) = IDloads(:,cycle,3); % Knee Flexion / Extension
            b(:,4) = IDloads(:,cycle,4); % Ankle Flexion / Extension
            b(:,5) = IDloads(:,cycle,5); % Ankle Inversion / Eversion

%             Comment plotting
%             % Plot moments
%             figure
%             subplot(1,5,1);
%             hold on;
%             plot(b(:,1),'LineWidth',1,'Color','k');
%             title ('Hip Flexion / Extension');
%             xlabel('% of Gait Cycle');
%             ylabel('Inter-segmental Moment (in N/m)');
%             subplot(1,5,2);
%             hold on;
%             plot(b(:,2),'LineWidth',1,'Color','k');
%             title ('Hip Adduction / Abduction');
%             xlabel('% of Gait Cycle');
%             subplot(1,5,3);
%             hold on;
%             plot(b(:,3),'LineWidth',1,'Color','k');
%             title ('Knee Flexion / Extension');
%             xlabel('% of Gait Cycle');
%             subplot(1,5,4);
%             hold on;
%             plot(b(:,4),'LineWidth',1,'Color','k');
%             title ('Ankle Flexion / Extension');
%             xlabel('% of Gait Cycle');
%             subplot(1,5,5);
%             hold on;
%             plot(b(:,5),'LineWidth',1,'Color','k');
%             title ('Ankle Inversion / Eversion');
%             xlabel('% of Gait Cycle');

            % Lever arms
            A = [permute(HipFEMA(:,cycle,:), [2,3,1]); ...
                permute(HipAAMA(:,cycle,:), [2,3,1]); ...
                permute(KneeMA(:,cycle,:), [2,3,1]); ...
                permute(AnkleMA(:,cycle,:), [2,3,1]); ...
                permute(SubtalarMA(:,cycle,:), [2,3,1])];

            % Muscle force
            f = permute(muscleForce(:,cycle,:), [3,2,1]);

            % Muscle passive force
            fpassive = permute(passiveF(:,cycle,:), [3,2,1]);

            % Muscle active force
            factive = f - fpassive;

            % Estimated moments from calibrated EMG
            for k = 1:101 % All frames (heel strike to heel strike)
%                 b2(:,:,k) = A(:,:,k)*f(:,:,k);
                b2(:,:,k) = A(:,:,k)*diag(cos(alpha))*f(:,:,k);
            end
            
%             Comment plotting
%             subplot(1,5,1); plot(squeeze(b2(1,1,:)),'LineWidth',1,'Color','b');
%             subplot(1,5,2); plot(squeeze(b2(2,1,:)),'LineWidth',1,'Color','b');
%             subplot(1,5,3); plot(squeeze(b2(3,1,:)),'LineWidth',1,'Color','b');
%             subplot(1,5,4); plot(squeeze(b2(4,1,:)),'LineWidth',1,'Color','b');
%             subplot(1,5,5); plot(squeeze(b2(5,1,:)),'LineWidth',1,'Color','b');

            % Activation
            activation = permute(a((cycle-1)*101+1:cycle*101,1:34),[2,3,1]);


            %% Muscle parameters

            % Ratio of slow-twitch fibers (Johnson et al. 1973)
            r =  [(0.535 + 0.633)/2; ... % addbrev (same as addmag)
                (0.535 + 0.633)/2; ... % addlong (same as addmag)
                (0.535 + 0.633)/2; ... % addmagDist
                (0.535 + 0.633)/2; ... % addmagIsch
                (0.535 + 0.633)/2; ... % addmagMid
                (0.535 + 0.633)/2; ... % addmagProx
                0.524; ... % glmax1
                0.524; ... % glmax2
                0.524; ... % glmax3
                0.524; ... % glmed1 (same as glmax)
                0.524; ... % glmed2 (same as glmax)
                0.524; ... % glmed3 (same as glmax)
                0.524; ... % glmin1 (same as glmax)
                0.524; ... % glmin2 (same as glmax)
                0.524; ... % glmin3 (same as glmax)
                0.492; ... % iliacus
                0.492; ... % psoas
                0.500; ... % tfl (taken from Uchida et al. (2016) model)
                0.669; ... % semimem (same as bf)
                0.669; ... % semiten  (same as bf)
                0.669; ... % bflh
                0.669; ... % bfsh
                (0.295 + 0.420 + 0.428)/3; ... % recfem
                (0.437 + 0.615)/2; ... % vasmed 
                (0.378 + 0.469)/2; ... % vaslat 
                (0.437 + 0.615 + 0.378 + 0.469)/4; ... % vasint (mean of vasmed & vaslat)
                (0.435 + 0.503 + 0.508)/3; ... % gaslat
                (0.435 + 0.503 + 0.508)/3; ... % gasmed
                (0.734 + 0.727)/2; ... % tibant
                (0.734 + 0.727)/2; ... % tibpost (same as tibant)
                0.625; ... % perbrev (same as perlong) 
                0.625; ... % perlong
                0.625; ... % pertert (same as perlong)
                (0.864 + 0.890)/2]; % soleus

            % Maximum isometric force (in N)
            % Peak isometric force values were calculated using information reported
            % in Handsfield et al. (2014)
            f0 = Fmax';

            % Muscle physiological cross sectional area (in m2)
            pcsa = f0/(61*1000); % Specific tension of 61 N/cm2

            % Muscle volume (in m3)
            volume = pcsa.*lmoOpt'; % PCSA * optimal muscle fibre length

            % Optimal length
            lmo = lmoOpt';
            % Slack length
            lts = ltsOpt';
            % Penation angle
            penation = alpha';

            % Muscle mass (in kg)
            % Unchida et al. (2016)
            mass = 1059.7*volume; % Density of mammalian muscle: 1059.7 kg/m3

            % Muscle-tendon length (in m)
            lmt = permute(Lmt(:,cycle,:), [3,2,1]);
            % Muscle-tendon velocity (in m/s)
            vmt = permute(Vmt(:,:,cycle), [2,3,1]);

            % Muscle force (from activation)
            % Maximal force when activation is 1
            % Minimal force when activation is 0 (passive)
            for j = 1:34
                for k = 1:101
%                     force_length_velocity(j,:,k) = factive(j,:,k) / ...
%                         (f0(j) * cos(alpha(j)) * activation(j,:,k));
%                     fmax(j,:,k) = f0(j) * cos(alpha(j)) * ...
%                         force_length_velocity(j,:,k) + ... % Activation set to 1
%                         fpassive(j,:,k);
                    force_length_velocity(j,:,k) = factive(j,:,k) / ...
                        (f0(j) * activation(j,:,k));
                    fmax(j,:,k) = f0(j) * ...
                        force_length_velocity(j,:,k) + ... % Activation set to 1
                        fpassive(j,:,k);
                    fmin(j,:,k) = fpassive(j,:,k); % Activation set to 0
                end
            end


            %% Cost functions

            % Sum of muscle forces
            J1 = sum(f,1); % At power 1
            J2 = sum(f.^2,1); % At power 2
            % J3 = sum(f.^3,1); % At power 3
            J4 = max(f);

            % Sum of activations
            a = (f-fmin)./fmax;
            J5 = sum(a,1); % At power 1
            J6 = sum(a.^2,1); % At power 2
            % J7 = sum(a.^3,1); % At power 3
            J8 = max(a);

            % Sum of muscle stresses
            stress = f./repmat(pcsa, [1 1 101]); % Muscle stress
            % (muscle stress and muscle force normalised by maximal isometric force are
            % similar criteria, i.e. just scaled by specific tension)
            J9 = sum(stress,1); % At power 1
            J10 = sum(stress.^2,1); % At power 2
            % J11 = sum(stress.^3,1); % At power 3
            J12 = max(stress);

            % Sum of muscle powers
            J13 = sum((f.*vmt).^2,1); % At power 2

            % Sum of muscle forces scaled by maximal muscle moments
            for j = 1:34 % Estimated maximal muscle moments
                for k = 1:101 % All frames (heel strike to heel strike)
                fmaxjk = zeros(34,1); % Initialisation
                fmaxjk(j,1) = fmax(j,1,k); % Instantaneous maximal force 
                    b0 = abs(A(:,:,k)*diag(cos(alpha))*fmaxjk); % Moment arm can be positive at one joint and negative at the other
                    M(j,:,k) = sum(b0)/sum(~b0 == 0); % Mean is case of biarticular muscle  
                end
            end
            J14 = sum((f./M).^2,1);

            % Others
            % Different normalisation
            fnormalisedf0 = f./repmat(f0, [1 1 101]); % Maximal isometric force
            factivenormalisedf0 = factive./repmat(f0, [1 1 101]); % Muscle force normalised by maximal isometric force
            stressactive = factive./repmat(pcsa, [1 1 101]);
            fnormalisedfmax = f./fmax; % Instantaneous maximal force 
            %
            J15 = sum(mass.*(0.5*factivenormalisedf0 + 0.5*stressactive.^2)); % Metabolic energy-related
            J16 = max((exp(3.48 + 0.169*r).* ...
                (100*fnormalisedf0).^(-0.5 - 0.036*r)) ...
                .^(-1)); % Minimum fatigue
            J17 = - real(sum(sqrt(1-(fnormalisedfmax).^2),1)); % Soft saturation
            % Real part even if fnormalisedfmax should theoritically not be > 1

            % Sum of joint stiffness
            % This cost function is based of force and force rate
            drdq = [permute(HipFE_drdq((cycle-1)*101+1:cycle*101,1:34),[3,2,1]); ...
                permute(HipAA_drdq((cycle-1)*101+1:cycle*101,1:34),[3,2,1]); ...
                permute(Knee_drdq((cycle-1)*101+1:cycle*101,1:34),[3,2,1]); ...
                permute(Ankle_drdq((cycle-1)*101+1:cycle*101,1:34),[3,2,1]); ...
                permute(Subtalar_drdq((cycle-1)*101+1:cycle*101,1:34),[3,2,1])];
            dfdl = permute(HipFE_dFdl(1:34,(cycle-1)*101+1:cycle*101),[1,3,2]);
            dfdt = dfdl.*vmt; % Force rate
            for k = 1:101
                K(:,:,k) = -(drdq(:,:,k) * f(:,:,k) - ...
                    (A(:,:,k)*diag(cos(alpha))).^2 * (dfdt(:,:,k)./vmt(:,:,k)));
            end
            J18 = -sum(K(1:4,:,:),1); % Maximised joint stiffness (no subtalar)            

            % Sum of joint reaction forces
            % Cannot be implemented with the available variables
            
%             Comment plotting
%             %% Comparison (scaled by maximum value)
%             figure
%             hold on
%             plot(squeeze(J1)/max(J1))
%             plot(squeeze(J2)/max(J2))
%             % plot(squeeze(J3)/max(J3))
%             plot(squeeze(J4)/max(J4))
%             plot(squeeze(J5)/max(J5))
%             plot(squeeze(J6)/max(J6))
%             % plot(squeeze(J7)/max(J7))
%             plot(squeeze(J8)/max(J8))
%             plot(squeeze(J9)/max(J9))
%             plot(squeeze(J10)/max(J10))
%             % plot(squeeze(J11)/max(J11))
%             plot(squeeze(J12)/max(J12))
%             % plot(squeeze(J13)/max(J13))
%             plot(squeeze(J14)/max(J14))
%             % plot(squeeze(J15)/max(J15))
%             plot(squeeze(J16)/max(J16))
%             plot(squeeze(J17)/34) % 34 muscles
%             plot(squeeze(J18)/max(-J18)) % Negative cost function
%             title ('Cost Functions');
%             xlabel('% of Gait Cycle');
%             ylabel('J/max(J)');
%             legend;

%             Comment plotting
%             %% Musculo-tendon force
%             figure
%             plot(permute(fnormalisedf0,[3,1,2]))
%             legend({'addbrev','addlong','addmagDist','addmagIsch','addmagMid', ...
%                 'addmagProx','glmax1','glmax2','glmax3','glmed1','glmed2','glmed3',...
%                 'glmin1','glmin2','glmin3','iliacus','psoas','semimem','semiten',...
%                 'bflh','bfsh','recfem','vasmed','vaslat','vasint','gaslat','gasmed',...
%                 'tibant','tibpost','perbrev','perlong','pertert','soleus','edl','fdl'})
%             title ('Musculo-tendon forces');
%             xlabel('% of Gait Cycle');
%             ylabel('f/f0 (peak isometric force)');
            
            %% Cost functions modified by me
            % Sum of musculo-tendon forces
%             J1 = sum(f,1) / n; % At power 1
%             J2 = sqrt(sum(f.^2,1) / n); % At power 2
%             J3 = (sum(f.^3,1) /n).^(1/3); % At power 3
%             J4 = max(f);
            J1 = zeros(1, 1, size(f, 3));
            J2 = zeros(1, 1, size(f, 3));
            J3 = zeros(1, 1, size(f, 3));
            J4 = zeros(1, 1, size(f, 3));            
            for k = 1 : size(f, 3)
                J1(:, k) = cost_function1(f(:, :, k));
                J2(:, k) = cost_function2(f(:, :, k));
                J3(:, k) = cost_function3(f(:, :, k));
                J4(:, k) = cost_function4(f(:, :, k));
            end
            % Sum of activations
%             a = (f-fmin)./fmax;
%             J5 = sum(a,1) / n; % At power 1
%             J6 = sqrt(sum(a.^2,1) / n); % At power 2
%             J7 = (sum(a.^3,1) / n).^(1/3); % At power 3
%             J8 = max(a);
            J5 = zeros(1, 1, size(f, 3));
            J6 = zeros(1, 1, size(f, 3));
            J7 = zeros(1, 1, size(f, 3));
            J8 = zeros(1, 1, size(f, 3));    
            for k = 1 : size(f, 3)
                J5(:, k) = cost_function5(f(:, :, k), fmin(:, :, k), fmax(:, :, k));
                J6(:, k) = cost_function6(f(:, :, k), fmin(:, :, k), fmax(:, :, k));
                J7(:, k) = cost_function7(f(:, :, k), fmin(:, :, k), fmax(:, :, k));
                J8(:, k) = cost_function8(f(:, :, k), fmin(:, :, k), fmax(:, :, k));
            end
            % Sum of muscle stresses
            stress = f./repmat(pcsa, [1 1 101]); % Muscle stress
            % (muscle stress and muscle force normalised by maximal isometric force are
            % similar criteria, i.e. just scaled by specific tension)
%             J9 = sum(stress,1) / n; % At power 1
%             J10 = sqrt(sum(stress.^2,1) / n); % At power 2
%             J11 = (sum(stress.^3,1) / n).^(1/3); % At power 3
%             J12 = max(stress);
            J9 = zeros(1, 1, size(f, 3));
            J10 = zeros(1, 1, size(f, 3));
            J11 = zeros(1, 1, size(f, 3));
            J12 = zeros(1, 1, size(f, 3));
            for k = 1 : size(f, 3)
                J9(:, k) = cost_function9(f(:, :, k), pcsa);
                J10(:, k) = cost_function10(f(:, :, k), pcsa);
                J11(:, k) = cost_function11(f(:, :, k), pcsa);
                J12(:, k) = cost_function12(f(:, :, k), pcsa);
            end
            % Sum of muscle powers
%             J13 = sqrt(sum((f.*vmt).^2,1) / n); % At power 2
            J13 = zeros(1, 1, size(f, 3));
            for k = 1 : size(f, 3)
                J13(:, k) = cost_function13(f(:, :, k), vmt(:, :, k));
            end
            
            % Sum of musculo-tendon forces scaled by maximal muscle moments
            for j = 1:34 % Estimated maximal muscle moments
                for k = 1:101 % All frames (heel strike to heel strike)
                fmaxjk = zeros(34,1); % Initialisation
                fmaxjk(j,1) = fmax(j,1,k); % Instantaneous maximal force 
                    b0 = abs(A(:,:,k)*diag(cos(alpha))*fmaxjk); % Moment arm can be positive at one joint and negative at the other
                    M(j,:,k) = sum(b0)/sum(~b0 == 0); % Mean is case of biarticular muscle  
                end
            end
%             J14 = sqrt(sum((f./M).^2,1) / n);
            J14 = zeros(1, 1, size(f, 3));
            for k = 1 : size(f, 3)
                J14(:, k) = cost_function14(f(:, :, k), M(:, :, k));
            end
            
            % Others
            % Different normalisation
            fnormalisedf0 = f./repmat(f0, [1 1 101]); % Maximal isometric force
            factivenormalisedf0 = factive./repmat(f0, [1 1 101]); % Muscle force normalised by maximal isometric force
            stressactive = factive./repmat(pcsa, [1 1 101]);
            fnormalisedfmax = f./fmax; % Instantaneous maximal force 
            %
%             J15 = sqrt( sum(mass.*(0.5*factivenormalisedf0 + 0.5*stressactive.^2)) / n); % Metabolic energy-related
%             J16 = max((exp(3.48 + 0.169*r).* ...
%                 (100*fnormalisedf0).^(-0.5 - 0.036*r)) ...
%                 .^(-1)); % Minimum fatigue
%             J17 = - real(sum(sqrt(1-(fnormalisedfmax).^2),1)) / n; % Soft saturation
            J15 = zeros(1, 1, size(f, 3));
            J16 = zeros(1, 1, size(f, 3));
            J17 = zeros(1, 1, size(f, 3));
            for k = 1 : size(f, 3)
                J15(:, k) = cost_function15(f(:, :, k), fpassive(:, :, k), f0, pcsa, mass);
                J16(:, k) = cost_function16(f(:, :, k), f0, r);
                J17(:, k) = cost_function17(f(:, :, k), fmax(:, :, k));
            end

            % Create cost function array
            J_arr = cat(1, J1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17);

            %% Store data
            % Time-variable data
            % Dimensionality same as force
            data.f(:, :, :, trial, speed-1, leg) = f;
            data.fmin(:, :, :, trial, speed-1, leg) = fmin;
            data.fmax(:, :, :, trial, speed-1, leg) = fmax;
            data.lmt(:, :, :, trial, speed-1, leg) = lmt;
            data.vmt(:, :, :, trial, speed-1, leg) = vmt;
            data.M(:, :, :, trial, speed-1, leg) = M;
            data.fpassive(:, :, :, trial, speed-1, leg) = fpassive;
            % Dimensionality same as moment
            data.A_orig(:, :, :, trial, speed-1, leg) = A;
            for k = 1 : 101
                data.A(:, :, k, trial, speed-1, leg) = A(:, :, k) * diag(cos(alpha));
            end
            data.b(:, :, :, trial, speed-1, leg) = b2;
            data.drdq(:, :, :, trial, speed-1, leg) = drdq;

            % Time-constant data
            % Dimensionality same as force
            data.pcsa(:, :, trial, speed-1, leg) = pcsa;
            data.f0(:, :, trial, speed-1, leg) = f0;
            data.mass(:, :, trial, speed-1, leg) = mass;
            data.r(:, :, trial, speed-1, leg) = r;
            data.lmo(:, :, trial, speed-1, leg) = lmo;
            data.lts(:, :, trial, speed-1, leg) = lts;
            data.penation(:, :, trial, speed-1, leg) = penation;
            
            
            % Cost function data            
            data.J_max = max(data.J_max, max(J_arr, [], 3).');
            data.J_min = min(data.J_min, min(J_arr, [], 3).');
            
%             % Computation of dfdl from activation
%             for j = 1:34
% 
%                 % Muscle parameters
%                 Params.Fmax_sym = f0(j);
%                 Params.lmo_sym = lmo(j);
%                 Params.lts_sym = lts(j);
%                 Params.alphaAngle = penation(j);
%                 Params.VmaxFactor_sym = 10;
%                 Params.Lmt_sym = squeeze(lmt(j,:,:));
%                 Params.Vmt_sym = squeeze(vmt(j,:,:));
%                 Params.a_sym = squeeze(a(j,:,:));
%                 % Computation of dfdl based on symbolic derivaties
%                 dfdl_sym(j,1,1:101) = permute(dfdlcal_v2(Params),[1,3,2]); 
% 
%             end
% 
%             %
%             % Stiffness estimated with the time derivative of force
%             for k = 1:101
%                 K_sym(:,:,k) = -(drdq(:,:,k) * f(:,:,k) - ...
%                     A(:,:,k).^2 * (dfdl_sym(:,:,k)));
%             end
            % Stiffness from the dataset
            K = [permute(Jstiffness_HipFE(1,(cycle-1)*101+1:cycle*101),[1,3,2]); ...
                permute(Jstiffness_HipAA(1,(cycle-1)*101+1:cycle*101),[1,3,2]); ...
                permute(Jstiffness_Knee(1,(cycle-1)*101+1:cycle*101),[1,3,2]); ...
                permute(Jstiffness_Ankle(1,(cycle-1)*101+1:cycle*101),[1,3,2]); ...
                permute(Jstiffness_Subtalar(1,(cycle-1)*101+1:cycle*101),[1,3,2])];
            
            data.K(:,:,:, trial, speed-1, leg) = K;
            fprintf('%d %d %d\n',trial,speed,leg);
            
            %% Foot-off
            % Peaks of extension (i.e. plantarlexion) moment
            [~,indminM] = min(b(:,4)); % Plantarlexion moment is negative
            % Just before the first positive value after the peak
            FootoffM = indminM + find(b(indminM + 1:end, 4)>0, 1, 'first') - 1;
            % Most of the time = peak of plantarflexion angle
            [~,indmina] = min(theta(indminM + 1:end, 4)); % Plantarlexion angle is negative
            Footoffa = indminM + indmina;
            
            % Store
            data.foot_off_index_m(:, :, trial, speed-1, leg) = FootoffM;
            data.foot_off_index_a(:, :, trial, speed-1, leg) = Footoffa;
            
            fprintf('Footoff: M = %d, A = %d\n', FootoffM, Footoffa);            
        end
    end
end

%% Save data

save_flag = input("Do you want to save data? [ANY KEY to save, ENTER to skip] : \n");
save_dir = "Optimization Model Data\Patient5.mat";
if save_flag
    save(save_dir, "data");    
end