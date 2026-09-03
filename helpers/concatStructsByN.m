function outStruct = concatStructsByN(structArray)
    % Validate input
    if ~iscell(structArray) || isempty(structArray)
        error('Input must be a non-empty cell array of structures.');
    end
    
    % Get common field names from the first structure
    fieldNames = fieldnames(structArray{1});
    
    % Initialize output structure
    outStruct = struct();
    for i = 1:numel(fieldNames)
        outStruct.(fieldNames{i}) = [];
    end
    
    % Loop through structures
    for i = 1:numel(structArray)
        S = structArray{i};
        
        % Basic consistency check
        if ~all(isfield(S, fieldNames))
            error('Structure %d does not have the expected fields.', i);
        end
        
        for j = 1:numel(fieldNames)
            field = fieldNames{j};
            val = S.(field);
            
            if isempty(outStruct.(field))
                outStruct.(field) = val;
            else
                % Check M-dimension match
                sz1 = size(outStruct.(field));
                sz2 = size(val);
                
                % Match along dimension 2 if both are 2D
                if numel(sz1) < 2, sz1(2) = 1; end
                if numel(sz2) < 2, sz2(2) = 1; end
                
                if sz1(2) ~= sz2(2)
                    error('Field "%s" has mismatched column size at entry %d.', field, i);
                end
                
                % Concatenate along N (rows)
                outStruct.(field) = cat(1, outStruct.(field), val);
            end
        end
    end
end
