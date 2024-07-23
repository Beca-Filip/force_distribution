function expanded_interval = expand_interval(interval, scale_factor)
%EXPAND_INTERVAL Expands an interval by a given scale factor.
%
% Syntax: expanded_interval = expand_interval(interval, scale_factor)
%
% Inputs:
%    interval - A 1x2 array representing the interval [start, end].
%    scale_factor - A scalar value by which to scale the length of the interval.
%
% Outputs:
%    expanded_interval - A 1x2 array representing the expanded interval.

    % Calculate the mean of the interval
    mean_interval = mean(interval);
    
    % Calculate the length of the interval
    interval_length = interval(2) - interval(1);
    
    % Calculate the expanded length
    expanded_length = scale_factor * interval_length;
    
    % Initialize the expanded interval with the same size and class as the input interval
    expanded_interval = zeros(size(interval), class(interval));
    
    % Calculate the expanded interval
    expanded_interval(1) = mean_interval - expanded_length / 2;
    expanded_interval(2) = mean_interval + expanded_length / 2;
end
