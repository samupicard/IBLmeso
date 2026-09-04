function argsOut = removeNameValues(argsIn, namesToRemove)
% removeNameValues
%
% Remove specified name/value pairs from a varargin-style cell array.

argsOut = {};
namesToRemove = string(namesToRemove);

i = 1;
while i <= numel(argsIn)
    if i < numel(argsIn) && (ischar(argsIn{i}) || isstring(argsIn{i}))
        nm = string(argsIn{i});
        if any(strcmpi(nm, namesToRemove))
            i = i + 2; % skip name and value
            continue
        end
    end
    argsOut{end+1} = argsIn{i}; %#ok<AGROW>
    i = i + 1;
end
end