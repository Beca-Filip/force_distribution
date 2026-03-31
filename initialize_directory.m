% Add casadi: find latest casadi directory for current OS in parent directory and add to path
parentDir = fileparts(mfilename('fullpath')); % directory containing this script
if isempty(parentDir)
    parentDir = pwd;
end
parentDir = fullfile(parentDir, '..'); % parent directory

% Determine OS identifier used in casadi folder names
if ispc
    osPattern = 'windows';
elseif ismac
    osPattern = 'apple-darwin';
elseif isunix
    osPattern = 'linux';
else
    error('Unsupported OS');
end

% List directories in parentDir matching casadi-* and containing osPattern
d = dir(parentDir);
isDir = [d.isdir];
dirNames = {d(isDir).name};
matches = {};
for k = 1:numel(dirNames)
    name = dirNames{k};
    if startsWith(name, 'casadi-','IgnoreCase',true) && contains(name, osPattern, 'IgnoreCase',true)
        matches{end+1} = fullfile(parentDir, name); %#ok<AGROW>
    end
end

if isempty(matches)
    warning('No casadi directory found for OS pattern "%s" in parent directory "%s".', osPattern, parentDir);
else
    % Choose the latest by datenum of folder (modification time) or by version-like name if available
    infos = cellfun(@(p) dir(p), matches, 'UniformOutput', false);
    modTimes = zeros(size(matches));
    for k = 1:numel(matches)
        info = dir(matches{k});
        if isempty(info)
            modTimes(k) = 0;
        else
            modTimes(k) = max([info.datenum]);
        end
    end
    [~, idx] = max(modTimes);
    casadi_path = matches{idx};
    disp(casadi_path)
    addpath(casadi_path);
end

% Add all paths in this directory
this_dir = '../force_distribution';
addpath(genpath(this_dir));

% Add tile figures
addpath('../TileFigures/');