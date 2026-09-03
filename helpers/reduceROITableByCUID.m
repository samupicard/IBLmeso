function Tout = reduceROITableByCUID(T)
% recucedROITableByCUID  Reduce aggregate ROI table to one row per non-empty cUID.
%
% Empty cUID rows are kept as independent rows and are not matched together.
% date, session and FOV columns are dropped.
% Numeric/logical fields are averaged within each cUID group.

if isempty(T)
    Tout = removevars(T, intersect({'date','session','FOV'}, T.Properties.VariableNames));
    return;
end

% Ensure cUID is string
T.cUID = string(T.cUID);

hasCUID = T.cUID ~= "";

% Split rows
T_chronic = T(hasCUID,:);
T_empty   = T(~hasCUID,:);

% Drop date/session
dropVars = intersect({'date','session','FOV'}, T.Properties.VariableNames);
T_chronic = removevars(T_chronic, dropVars);
T_empty   = removevars(T_empty, dropVars);

% Empty cUIDs are kept as-is
if isempty(T_chronic)
    Tout = T_empty;
    return;
end

% Variables to average
varNames = T_chronic.Properties.VariableNames;
keyVars = {'cUID'};
nonMeanVars = {'subject','FOV','cUID'};

meanVars = setdiff(varNames, nonMeanVars, 'stable');

% Keep only numeric/logical variables for averaging
isMeanVar = false(size(meanVars));
for i = 1:numel(meanVars)
    x = T_chronic.(meanVars{i});
    isMeanVar(i) = isnumeric(x) || islogical(x);
end
meanVars = meanVars(isMeanVar);

% Group by cUID
[G, cUIDs] = findgroups(T_chronic.cUID);

nG = numel(cUIDs);

Tout_chronic = table;
Tout_chronic.cUID = cUIDs;

% subject: keep first value per cUID
if ismember('subject', varNames)
    Tout_chronic.subject = splitapply(@(x) x(1), T_chronic.subject, G);
end

% FOV: keep first value per cUID
if ismember('FOV', varNames)
    Tout_chronic.FOV = splitapply(@(x) x(1), T_chronic.FOV, G);
end

% Means for numeric/logical variables
for i = 1:numel(meanVars)
    v = meanVars{i};
    Tout_chronic.(v) = splitapply(@(x) mean(x, 'omitnan'), T_chronic.(v), G);
end

% Match original variable order after dropping date/session
outVars = setdiff(T.Properties.VariableNames, dropVars, 'stable');
outVars = intersect(outVars, Tout_chronic.Properties.VariableNames, 'stable');
Tout_chronic = Tout_chronic(:, outVars);

% Combine chronic-reduced rows with empty-cUID rows
Tout = [Tout_chronic; T_empty(:, Tout_chronic.Properties.VariableNames)];

end