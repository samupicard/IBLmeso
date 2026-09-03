function [orderedLabels,orderIndex,groupIndex] = ...
    orderAllenRegionAcronyms(regionLabels)
%ORDERALLENREGIONACRONYMS Apply a preferred Allen-region plotting order.
%
% Ordering:
%
%   1. VISp*
%   2. Other VIS* regions, excluding VISa* and VISrl*
%   3. VISa* and VISrl*
%   4. RSP*
%   5. AUD*
%   6. TE*
%   7. SS*
%   8. MO*
%   9. All remaining regions
%
% Within each group, acronyms are sorted alphabetically.
%
% INPUT
%   regionLabels
%       Vector of Allen atlas acronyms.
%
% OUTPUTS
%   orderedLabels
%       Region labels in preferred order.
%
%   orderIndex
%       Indices such that:
%
%           orderedLabels = regionLabels(orderIndex)
%
%   groupIndex
%       Ordering group for each element of orderedLabels.

arguments
    regionLabels
end

inputWasRow = isrow(regionLabels);

labels = string(regionLabels);
labels = labels(:);

nRegions = numel(labels);

groupForOriginal = 9*ones(nRegions,1);

%% Define wildcard-style groups

isVISp = startsWith( ...
    labels, ...
    "VISp", ...
    "IgnoreCase",true);

isVISa = ~cellfun( ...
    @isempty, ...
    regexpi( ...
        cellstr(labels), ...
        '^VISa(?:\d|$)', ...
        'once'));

isVISrl = startsWith( ...
    labels, ...
    "VISrl", ...
    "IgnoreCase",true);

isVIS = startsWith( ...
    labels, ...
    "VIS", ...
    "IgnoreCase",true);

%% Assign groups

groupForOriginal(isVISp) = 1;

groupForOriginal( ...
    isVIS & ...
    ~isVISp & ...
    ~isVISa & ...
    ~isVISrl) = 2;

groupForOriginal(isVISa | isVISrl) = 3;

groupForOriginal( ...
    startsWith(labels,"RSP","IgnoreCase",true)) = 4;

groupForOriginal( ...
    startsWith(labels,"AUD","IgnoreCase",true)) = 5;

groupForOriginal( ...
    startsWith(labels,"TE","IgnoreCase",true)) = 6;

groupForOriginal( ...
    startsWith(labels,"SS","IgnoreCase",true)) = 7;

groupForOriginal( ...
    startsWith(labels,"MO","IgnoreCase",true)) = 8;

%% Sort by group, then alphabetically within group

sortTable = table( ...
    groupForOriginal, ...
    upper(labels), ...
    (1:nRegions).');

sortTable.Properties.VariableNames = { ...
    'group', ...
    'label', ...
    'originalIndex'};

sortTable = sortrows( ...
    sortTable, ...
    {'group','label'});

orderIndex = sortTable.originalIndex;
groupIndex = groupForOriginal(orderIndex);
orderedLabels = labels(orderIndex);

%% Preserve input orientation

if inputWasRow
    orderedLabels = orderedLabels.';
    orderIndex = orderIndex.';
    groupIndex = groupIndex.';
end

end