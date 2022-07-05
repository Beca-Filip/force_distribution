function data = compute_stiffness_wrap(data, sample_list, trial_list, speed_list, leg_list)

% Number of optimizations to do
nsamples = length(sample_list);
ntrials = length(trial_list);
nspeeds = length(speed_list);
nlegs = length(leg_list);

% % Create output structure
% Kout = zeros([size(data.f,[1 2]), nsamples, ntrials, nspeeds, nlegs]);


% Calculate activation
a = (data.f - data.fmin) ./ (data.fmax - data.fmin);

% Counters
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
            
%             % Reset sample counter
%             cntsamples = 1;
%             % Loop over samples
%             for k = sample_list

                
                % Predict stiffness
                % Computation of dfdl from activation
                for j = 1:35

                    % Muscle parameters
                    Params.Fmax_sym = data.f0(j, :, trial, speed, leg);
                    Params.lmo_sym = data.lmo(j, :, trial, speed, leg);
                    Params.lts_sym = data.lts(j, :, trial, speed, leg);
                    Params.alphaAngle = data.penation(j, :, trial, speed, leg);
                    Params.VmaxFactor_sym = 10;
                    Params.Lmt_sym = squeeze(data.lmt(j,:,sample_list,trial,speed,leg));
                    Params.Vmt_sym = squeeze(data.vmt(j,:,sample_list,trial,speed,leg));
                    Params.a_sym = squeeze(a(j,:,sample_list,trial,speed,leg));
                    % Computation of dfdl based on symbolic derivaties
                    dfdl_sym(j,1,1:nsamples) = permute(dfdlcal_v2(Params),[1,3,2]); 
                end
                
                % Stiffness estimated with the length derivative of force
                for k = 1:nsamples
                    K_sym(:,:,k) = -(data.drdq(:,:,sample_list(k),trial,speed,leg) * data.f(:,:,sample_list(k),trial,speed,leg) - ...
                        data.A(:,:,sample_list(k),trial,speed,leg).^2 * (dfdl_sym(:,:,k)));
                end
                
                % Store stiffness
                data.K(:,:,sample_list,trial,speed,leg)=K_sym;

%                 % Augment sample counter
%                 cntsamples = cntsamples + 1;
%             end
            
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