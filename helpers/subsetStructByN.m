function outStruct = subsetStructByN(inStruct, idx)
    % Validate input
    if islogical(idx)
        isLogical = true;
        nExpected = numel(idx);  % Used to match against field sizes
    elseif isnumeric(idx)
        isLogical = false;
    else
        error('Index must be a logical array or a vector of integers.');
    end

    fieldNames = fieldnames(inStruct);
    outStruct = inStruct;

    for i = 1:numel(fieldNames)
        field = fieldNames{i};
        val = inStruct.(field);
        
        if isempty(val)
            continue;
        end
        
        sz = size(val);
        
        %assume first dimension of first field is size N
        if i==1 && ~isLogical
            nExpected = sz(1);
        end
        
        % Check if first dimension matches N
        if sz(1) == nExpected
            outStruct.(field) = val(idx, :);
        
        % Check if second dimension matches N (e.g., row vector fields)
        elseif sz(2) == nExpected
            outStruct.(field) = val(:, idx);
        
        else
            % Field does not match size N, leave unchanged
        end
    end
end
