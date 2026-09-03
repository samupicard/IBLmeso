function [avgResps, typeVals] = getAvgResps(resps, trialsT, typeFld, typeVals, avgType)

% getAvgResps computes average responses (PETH) per trial type
%
%   inputs
%     resps     - nTrials x nTimepoints, matrix of neuron's responses
%     trialsT   - nTrials x nCols, trials table with trial types and event latencies
%     typeFld   - string. Trial type to group by (e.g. contrastDiff)
%     typeVals  - (optional) values of typeFld to select
%     avgType   - (optional) 'mean' (default) or 'median'
%
%   outputs
%     avgResps  - nTypeVals x nTimepoints matrix. Average response for each
%                 trial type specified in typeVals.
%     typeVals  - vector of trial-type values corresponding to rows of
%                 avgResps.
%
%   Samuel Picard (Jan 2026)

% --- defaults ---
if nargin < 4 || isempty(typeVals)
    vals = trialsT{:,typeFld};
    typeVals = unique(vals);
else
    vals = trialsT{:,typeFld};
end

if nargin < 5 || isempty(avgType)
    avgType = 'mean';
end

% --- select averaging function ---
switch lower(avgType)
    case 'mean'
        avgFun = @(x) nanmean(x,1);
    case 'median'
        avgFun = @(x) nanmedian(x,1);
    otherwise
        error('avgType must be ''mean'' or ''median''.');
end

% --- indices per trial type ---
iV = arrayfun(@(v) find(vals == v), typeVals, 'UniformOutput', false);

% --- compute averages ---
avgResps = nan(length(typeVals), size(resps,2));
for iVal = 1:length(iV)
    avgResps(iVal,:) = avgFun(resps(iV{iVal},:));
end
