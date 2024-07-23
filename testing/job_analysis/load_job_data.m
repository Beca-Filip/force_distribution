function trials = load_job_data(dirname)
    % Get a list of all files in the specified directory
    files = dir(fullfile(dirname, '*.mat'));
    
    trials = {};
    % Loop over each file in the directory
    for i = 1:length(files)
        % Get the full path to the current file
        filepath = fullfile(dirname, files(i).name);
        
        % Extract the integers i, j, and k from the filename
        tokens = regexp(files(i).name, '(\d+)-(\d+)-(\d+)', 'tokens');
        
        if ~isempty(tokens)
            phase_sel = str2double(tokens{1}{1});
            speed_sel = str2double(tokens{1}{2});
            leg_sel = str2double(tokens{1}{3});
            
            % Import data from the current file
            trials{leg_sel, speed_sel, phase_sel} = importdata(filepath);
            
            % Process the data as needed
            % (e.g., store it in a variable, perform analysis, etc.)
            fprintf('Processing trial: %s with phase=%d, speed=%d, leg=%d\n', files(i).name, phase_sel, speed_sel, leg_sel);
            % Example: display the imported data
            % disp(trial);
        else
            warning('Filename %s does not match the expected pattern.', files(i).name);
        end
    end
end
