function commonLabels = findCommonLabels(labels_bySession, session_selection)

    % Validate input
    if ~iscell(labels_bySession) || isempty(labels_bySession)
        error('Input must be a non-empty cell array.');
    end
    
    N = numel(labels_bySession);
    
    % Default: use all cells
    if nargin < 2 || isempty(session_selection)
        session_selection = true(1, N);
    end
    
    % Ensure selection is logical
    if isnumeric(session_selection)
        tmp = false(1, N);
        tmp(session_selection) = true;
        session_selection = tmp;
    elseif ~islogical(session_selection) || numel(session_selection) ~= N
        error('Selection must be a logical of the same length as the cell array, or a numeric vector.');
    end
    
    selectedCells = labels_bySession(session_selection);
    
    % Start with the unique elements of the first selected cell
    commonLabels = unique(selectedCells{1});
    
    % Intersect with the rest
    for i = 2:length(selectedCells)
        commonLabels = intersect(commonLabels, unique(selectedCells{i}));
        if isempty(commonLabels)
            break;
        end
    end
end
