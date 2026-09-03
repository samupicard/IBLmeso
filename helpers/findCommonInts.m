function commonVals = findCommonInts(cellArray)
    % Validate input
    if ~iscell(cellArray) || isempty(cellArray)
        error('Input must be a non-empty cell array.');
    end
    
    % Start with the unique elements of the first cell
    commonVals = unique(cellArray{1});
    
    % Intersect with each subsequent cell's contents
    for i = 2:length(cellArray)
        commonVals = intersect(commonVals, unique(cellArray{i}));
        
        % Early exit: if nothing left in common
        if isempty(commonVals)
            break;
        end
    end
end
