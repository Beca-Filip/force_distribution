close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Filtered_Patient5.mat';
load(data_dir);

% Modify normalization
data.J_min(:) = 0;
data.J_max = data.J_max ./ 1e3;

% Exclude cf
cf_exclude = [16, 17];

% Create model
[model, vars] = form_casadi_model_normalized_s5(cf_exclude);

% Choose solver with options
sol_opt= struct;
sol_opt.ipopt.print_level = 0;
sol_opt.print_time =0;
sol_opt.verbose = 0;
sol_opt.ipopt.sb ='yes';
sol_opt.ipopt.check_derivatives_for_naninf = 'yes';
sol_opt.regularity_check = true;
model.solver('ipopt', sol_opt);

%% Perform IOC on these (To compare with RMSE found in test_workspace_do_1)
nspeeds = 5;
nlegs = 2;
nphases = 2;
trial_list = 1:10;

% Sample list
sample_cell_list = cell(nphases, length(trial_list), nspeeds, nlegs);

for speed = 1 : 5
    for leg = 1 : 2
        for phase = 1 : 2
            
            % If it's the stance phase
            if phase == 1
                
                % For each trial take every fourth sample until
                % the foot is taken off the ground
                for trial = trial_list
                    % takeoff sample
                    sample_cell_list{phase, trial, speed, leg} = 1 : 4 : data.foot_off_index_m(:,:,trial,speed,leg);
                end
                
            % If it's the swing phase
            elseif phase == 2
                
                % For each trial take every fourth sample starting from the
                % sample when the foot is lifted off
                for trial = trial_list
                    % takeoff sample
                    sample_cell_list{phase, trial, speed, leg} = data.foot_off_index_m(:,:,trial,speed,leg) : 4 : 101;
                end
            end
        end
    end
end


for speed = 1 : 5
    for leg = 1 : 2
        for phase = 1 : 2
            
            fprintf("Working on: speed %d, leg %d, phase %d.\n", speed, leg, phase);
            % Get number of grid points
            Ngrid = 3060;
            
            % Calculate RMSE
            tic
            [E, alpha] = IO_grid_search_normalized_variable_phase(Ngrid, data, vars, model, sample_cell_list, trial_list, speed, leg);
            toc
            
            % Saving 
            prepresuffix = sprintf('Ngrid_%d', Ngrid);
            presuffix = sprintf('leg_%d_speed_%d_phase_%d_', leg, speed, phase);
            suffix = datetimestr;
            
            save(sprintf('..\\..\\opti_results\\job_patient5\\diffnorm_gs_%s_%s_%s.mat', prepresuffix, presuffix, suffix));
        end
    end
end