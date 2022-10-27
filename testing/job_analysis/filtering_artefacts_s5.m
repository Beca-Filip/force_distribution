close all;
clear all;
clc;

% Data directory and loading
data_dir = '..\..\Optimization Model Data\Patient5.mat';
load(data_dir);

% Define additional field: optimization initialization
% Should contain initialization values for the optimization model
data.finit_feas = zeros(size(data.f));

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

TAB = table();

sample_list = 0:100;

save_flag = 0;

% Plotting options
opts.xlabel = @(ii) {ternary_operator(ii >= 29, "$t$ [\%]", ""), 'interpreter', 'latex'};
opts.ylabel = @(ii) {sprintf("$f_{%d}$ [N]", ii), 'interpreter', 'latex'};

% Do every leg, speed and phase
for leg = [1, 2]
for speed = [1 : 5]
for trial = 1:10
    
    % Prepare gridding
    gridShape = [6, 6];
    fig = figure;
    
    % Figure out the title
    stitle = sprintf("Leg: %d, Speed: %d, Trial: %d", leg, speed, trial);
    if any(any(any( ( data.fmax(:,:,:,trial,speed,leg) - data.fmin(:,:,:,trial,speed,leg) ) <= 0)))
        stitle = strcat(stitle, "; Has fmin<=fmax");
    end
    if any(any(any( ( data.fmax(:,:,:,trial,speed,leg) - data.f(:,:,:,trial,speed,leg) ) < 0)))
        stitle = strcat(stitle, "; Has fmax < f");
    end
    if any(any(any( ( data.f(:,:,:,trial,speed,leg) - data.fmin(:,:,:,trial,speed,leg) ) < 0)))
        stitle = strcat(stitle, "; Has f < fmin");
    end
    suptitle(stitle);
    
    % Median filter the fmax <= fmin
    filtered_flag = 1;
    for k = 1 : size(data.f, 1)
        if any( ( data.fmax(k,:,:,trial,speed,leg) - data.fmin(k,:,:,trial,speed,leg) ) <= 0)
            zeroIndices = find(( data.fmax(k,:,:,trial,speed,leg) - data.fmin(k,:,:,trial,speed,leg) ) <= 0);
            for filtIndex = zeroIndices
                % median filter
                filterWidth = 20;
                filtBegin = max(1, filtIndex - filterWidth/2);
                filtEnd = min(size(data.fmax, 3), filtIndex + filterWidth/2);
                data.fmax(k, :, filtIndex, trial,speed,leg) = median(data.fmax(k, :, filtBegin : filtEnd, trial,speed,leg));
            end
        end
    end

    set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
    plot_vector_quantities_opts_shape(sample_list, data.f(:,:,:,trial,speed,leg), [], opts, gridShape, 'LineWidth', 2, 'Color', [0, 0, 1], 'DisplayName', 'Forces');
    plot_vector_quantities_opts_shape(sample_list, data.fmin(:,:,:,trial,speed,leg), [], [], gridShape, 'LineWidth', 2, 'Color', [0.3, 0.3, 0.3], 'DisplayName', 'fmin');
    plot_vector_quantities_opts_shape(sample_list, data.fmax(:,:,:,trial,speed,leg), [], [], gridShape, 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'fmax');
    
    
    % Plot the fmin >= fmax
    for k = 1 : size(data.f, 1)
        if any( ( data.fmax(k,:,:,trial,speed,leg) - data.fmin(k,:,:,trial,speed,leg) ) <= 0)
            subplot(gridShape(1), gridShape(2), k)
            hold on;
            zeroIndices = find(( data.fmax(k,:,:,trial,speed,leg) - data.fmin(k,:,:,trial,speed,leg) ) <= 0);
            scatter(zeroIndices-1, data.fmin(k,:,zeroIndices,trial,speed,leg), 20, 'r', 'DisplayName', '$f_{\rm max} <= f_{\rm min}$');
        end
    end
    
    % Plot the fmax <= f
    for k = 1 : size(data.f, 1)
        if any( ( data.fmax(k,:,:,trial,speed,leg) - data.f(k,:,:,trial,speed,leg) ) < 0)
            subplot(gridShape(1), gridShape(2), k)
            hold on;
            zeroIndices = find(( data.fmax(k,:,:,trial,speed,leg) - data.f(k,:,:,trial,speed,leg) ) < 0);
            scatter(zeroIndices-1, data.fmax(k,:,zeroIndices,trial,speed,leg), 20, 'g', 'DisplayName', '$f_{\rm max} <= f$');
        end
    end
    
    % Plot the fmax <= f
    for k = 1 : size(data.f, 1)
        if any( ( data.f(k,:,:,trial,speed,leg) - data.fmin(k,:,:,trial,speed,leg) ) < 0)
            subplot(gridShape(1), gridShape(2), k)
            hold on;
            zeroIndices = find(( data.fmax(k,:,:,trial,speed,leg) - data.f(k,:,:,trial,speed,leg) ) < 0);
            scatter(zeroIndices-1, data.fmax(k,:,zeroIndices,trial,speed,leg), 20, 'c', 'DisplayName', '$f_{\rm max} <= f$');
        end
    end
    
    % If f >= fmax, or f <= fmin, meaning if the measured f is infeasible
    % calculate the closest feasible f and store it into finit_feas
    % Otherwise, if it is feasible just store f into finit_feas
    for k = 1 : size(data.f, 3)
        % If current f infeasible
       if any( ( data.fmax(:,:,k,trial,speed,leg) - data.f(:,:,k,trial,speed,leg) ) < 0) || ...
          any( ( data.f(:,:,k,trial,speed,leg) - data.fmin(:,:,k,trial,speed,leg) ) < 0)
        
           % Calculate the closest feasiblepoint to the current force
           ineqSlack = .05 * (data.fmax(:,:,k,trial,speed,leg) - data.fmin(:,:,k,trial,speed,leg));
           data.finit_feas(:, :, k, trial, speed, leg) = ...
           lsqlin(eye(size(data.f, 1)),data.f(:,:,k,trial,speed,leg),[],[], ...
           data.A(:,:,k,trial,speed,leg),data.b(:,:,k,trial,speed,leg),....
           data.fmin(:,:,k,trial,speed,leg)+ineqSlack,data.fmax(:,:,k,trial,speed,leg)-ineqSlack,[],...
           optimoptions(@lsqlin, 'Display', 'final', 'Algorithm', 'interior-point', 'ConstraintTolerance', 1e-10));
       else
           % Just store the measured f
           data.finit_feas(:, :, k, trial, speed, leg) = data.f(:, :, k, trial, speed, leg);
       end
    end
    
    % Plot feasible points
    plot_vector_quantities_opts_shape(sample_list, data.finit_feas(:,:,:,trial,speed,leg), [], [], gridShape, '--', 'LineWidth', 2, 'Color', [1, 0.2, 0.2], 'DisplayName', 'ffeas');
    
    if save_flag
        if ~filtered_flag
            saveas(gcf, sprintf('temp/trj_leg_%d_speed_%d_trial_%d.png', leg, speed, trial));
        else
            saveas(gcf, sprintf('temp/filt_trj_leg_%d_speed_%d_trial_%d.png', leg, speed, trial));
        end
    end
end    
end
end

% Save
if filtered_flag
    save_data_dir = '..\..\Optimization Model Data\Filtered_Patient5.mat';
    save(save_data_dir, 'data');
end