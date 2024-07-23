function propArray = extract_properties_from_structs(cellArray, propName)
    % Get the number of dimensions and size of the input cell array
    numDims = ndims(cellArray);
    arraySize = size(cellArray);

    % Get the size of the specified property in the first struct
    propSize = size(cellArray{1}.(propName));
    propDims = ndims(cellArray{1}.(propName));
    combinedSize = [arraySize, propSize];
    
    % Initialize the output array with NaN values
    propArray = nan(combinedSize);

    % Loop through each element in the cell array
    for idx = 1:numel(cellArray)
        % Get the indices of the current element in the cell array
        [propArrayIdx{1:numDims}] = ind2sub(arraySize, idx);
        

        % Construct the indexing expression for propArray
        for jj = 1 : propDims
            propArrayIdx{jj + numDims} = 1 : propSize(jj);
        end

        % Extract and assign the specified property from the current struct
        propArray(propArrayIdx{:}) = cellArray{idx}.(propName);
    end
end
