
function T_clean = removeEmptyRows(T,ignoreFields)

if nargin<2
    ignoreFields = {'annot','tstat_pseudo','tstat_empirical'};
end

vars      = T.Properties.VariableNames;
checkVars = setdiff(vars, ignoreFields);

nRows = height(T);
nVars = numel(checkVars);

% Preallocate logical matrix: rows x varsToCheck
isEmptyMat = false(nRows, nVars);

for j = 1:nVars
    col = T.(checkVars{j});
    if iscell(col)
        % cell columns: empty if the cell contents are empty
        isEmptyMat(:, j) = cellfun(@isempty, col);
        
    elseif isnumeric(col) || islogical(col)
        % numeric / logical: treat NaN as empty
        if size(col,1)>1 && size(col,2)>1
            isEmptyMat(:, j) = isnan(col(:,1));
        else
            isEmptyMat(:, j) = isnan(col);
        end
        
    elseif isstring(col)
        % string array: empty string
        isEmptyMat(:, j) = (col == "");
        
    elseif isdatetime(col) || isduration(col)
        % datetime/duration: NaT / missing
        isEmptyMat(:, j) = isnat(col);
        
    else
        % fallback: per-element isempty (for weird types)
        isEmptyMat(:, j) = arrayfun(@isempty, col);
    end
end

% Row is "empty" if *all* checked vars are empty
rowsToDrop = all(isEmptyMat, 2);

% Keep only rows where at least one non-annot field is non-empty
T_clean = T(~rowsToDrop, :);