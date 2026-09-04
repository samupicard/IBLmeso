function S = restrict_session_to_uids(S, uidList)
% restrict_session_to_uids
%
% Keep only rows whose chronicUID is in uidList, and reorder rows to match
% the order of uidList.
%
% Assumes S is a session struct from load_sessionPETH_all_data with rows
% already aligned across all peth types.
%
% Required fields:
%   S.chronicUID         [nROI x 1 string]
%
% Optional per-type fields that will be restricted if present:
%   S.Diff_even          {1 x nTypes}, each [nROI x nTime]
%   S.Diff_odd           {1 x nTypes}, each [nROI x nTime]
%   S.MeanPos_even       {1 x nTypes}, each [nROI x nTime]
%   S.MeanNeg_even       {1 x nTypes}, each [nROI x nTime]
%   S.MeanPos_odd        {1 x nTypes}, each [nROI x nTime]
%   S.MeanNeg_odd        {1 x nTypes}, each [nROI x nTime]
%   S.MI_even            {1 x nTypes}, each [nROI x nTime]
%   S.MI_odd             {1 x nTypes}, each [nROI x nTime]
%   S.sortMetric         {1 x nTypes}, each [nROI x 1]
%
% Optional row-wise fields that will also be restricted if present:
%   S.globalROI
%   S.cellScore
%   S.brainIds
%   S.fovLabel
%   S.isChronic
%   S.isResponsive
%   S.roiScale
%   S.sharedROI
%   S.mlapdv
%
% Example:
%   commonUID = intersect(S1.chronicUID, S2.chronicUID, 'stable');
%   commonUID = commonUID(strlength(commonUID) > 0);
%   S1 = restrict_session_to_uids(S1, commonUID);
%   S2 = restrict_session_to_uids(S2, commonUID);

uidList = string(uidList(:));
uidList = uidList(strlength(strtrim(uidList)) > 0);

if ~isfield(S, 'chronicUID')
    error('Session struct must contain field S.chronicUID');
end

S.chronicUID = string(S.chronicUID(:));

if isempty(uidList)
    error('uidList is empty after removing blank entries.');
end

% map requested UID order onto session rows
[tf, loc] = ismember(uidList, S.chronicUID);
rowIdx = loc(tf);

if isempty(rowIdx)
    error('None of the requested UIDs were found in S.chronicUID.');
end

% sanity: duplicated chronicUIDs inside session would make mapping ambiguous
u = unique(S.chronicUID(strlength(strtrim(S.chronicUID)) > 0));
if numel(u) ~= numel(S.chronicUID(strlength(strtrim(S.chronicUID)) > 0))
    warning('S.chronicUID contains duplicates; restrict_session_to_uids uses first-match behavior via ismember.');
end

% reorder row-wise metadata
rowFields = { ...
    'chronicUID', ...
    'globalROI', ...
    'cellScore', ...
    'brainIds', ...
    'fovLabel', ...
    'isChronic', ...
    'isResponsive', ...
    'roiScale', ...
    'sharedROI'};

for iF = 1:numel(rowFields)
    fn = rowFields{iF};
    if isfield(S, fn) && ~isempty(S.(fn))
        try
            S.(fn) = S.(fn)(rowIdx, :);
        catch
            % ignore fields that are not simple row vectors
        end
    end
end

% mlapdv is row-wise but 2D
if isfield(S, 'mlapdv') && ~isempty(S.mlapdv)
    S.mlapdv = S.mlapdv(rowIdx, :);
end

% reorder per-type matrices / vectors
typeFields = { ...
    'Diff_even', ...
    'Diff_odd', ...
    'MeanPos_even', ...
    'MeanNeg_even', ...
    'MeanPos_odd', ...
    'MeanNeg_odd', ...
    'MI_even', ...
    'MI_odd', ...
    'sortMetric'};

for iF = 1:numel(typeFields)
    fn = typeFields{iF};
    if isfield(S, fn) && ~isempty(S.(fn))
        for s = 1:numel(S.(fn))
            if ~isempty(S.(fn){s})
                S.(fn){s} = S.(fn){s}(rowIdx, :);
            end
        end
    end
end

% reset sharedROI to current row IDs if desired; here keep semantic UID order
if isfield(S, 'sharedROI')
    S.sharedROI = (1:numel(rowIdx))';
end
end