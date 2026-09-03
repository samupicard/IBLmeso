function outStruct = concatStructsByT(structArray)
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
            
            if j==1
                % we can assume first field has size [N x T]
                sz = size(val);
            end
            
            if isempty(outStruct.(field))
                outStruct.(field) = val;
            else
                % Check N-dimension match
                sz1 = size(outStruct.(field));
                sz2 = size(val);
                
                % Match along dimension 2 if both are 2D
                if numel(sz1) < 2, sz1(2) = 1; end
                if numel(sz2) < 2, sz2(2) = 1; end
                
                if sz1(2) > 1 && sz1(1) ~= sz2(1)
                    error('Field "%s" has mismatched row size at entry %d.', field, i);
                end
                
                % Concatenate along T if other dimension is matching,
                % otherwise create a new dimension
                if ~ismember(sz2(1), [sz(2),1]) && ~ismember(sz2(2), [sz(2),1])
                    outStruct.(field) = cat(3, outStruct.(field), val);
                elseif sz1(1) == sz2(1)
                    outStruct.(field) = cat(2, outStruct.(field), val);
                else
                    outStruct.(field) = cat(1, outStruct.(field), val);
                end
            end
        end
    end
end
