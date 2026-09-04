function S = restrict_session_to_types(S, keepTypes)
% restrict_session_to_types
%
% Keep only the requested PETH types from a session struct S, preserving
% the order given in keepTypes.
%
% Expected fields in S (if present) that are indexed by peth type:
%   pethTypes
%   T_byType
%   Diff_even
%   Diff_odd
%   MeanPos_even
%   MeanNeg_even
%   MeanPos_odd
%   MeanNeg_odd
%   MI_even
%   MI_odd
%   Diff_plot
%   sortMetric
%   sortIdxByType
%   globalROI_byType
%   fov_byType
%
% Example:
%   commonTypes = intersect(S1.pethTypes, S2.pethTypes, 'stable');
%   S1 = restrict_session_to_types(S1, commonTypes);
%   S2 = restrict_session_to_types(S2, commonTypes);

keepTypes = cellstr(string(keepTypes));

if ~isfield(S, 'pethTypes')
    error('Session struct must contain field S.pethTypes');
end

pethTypes = cellstr(string(S.pethTypes));
[tf, loc] = ismember(keepTypes, pethTypes);
loc = loc(tf);
keepTypes = keepTypes(tf);

if isempty(loc)
    error('None of the requested types were found in S.pethTypes.');
end

% preserve requested order
S.pethTypes = pethTypes(loc);

typeFields = { ...
    'T_byType', ...
    'Diff_even', ...
    'Diff_odd', ...
    'MeanPos_even', ...
    'MeanNeg_even', ...
    'MeanPos_odd', ...
    'MeanNeg_odd', ...
    'MI_even', ...
    'MI_odd', ...
    'Diff_plot', ...
    'sortMetric', ...
    'sortIdxByType', ...
    'globalROI_byType', ...
    'fov_byType'};

for iF = 1:numel(typeFields)
    fn = typeFields{iF};
    if isfield(S, fn) && ~isempty(S.(fn))
        S.(fn) = S.(fn)(loc);
    end
end
end