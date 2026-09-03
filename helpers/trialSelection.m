function trialCnd = trialSelection(trialsT,trialTypeField,trialTypeVals,trialTypeFilter)
%
% trialSelection makes a selection of trials from a IBL trials table
%
% inputs 
%   trialsT is the trials table (must be a table (not a struct) and must contain the trialTypeField)
%   trialTypeField defines trial type to look at e.g. 'contrastDiff'
%
% optional inputs
%   trialTypeVals is conditions to select e.g. [0.2, 0.8]
%   trialTypeFilter is a string; filters the trials according to a logical statement e.g. 'probabilityLeft<0.5'

% outputs
%   trialCnd logical array [nConditions, nTrials]
%
% Samuel Picard (Oct 2023)

if nargin<3
    trialTypeVals = '';
    trialTypeFilter = '';
elseif nargin<4
    trialTypeFilter = '';
end


% if a filter was defined, make this selection
if ~isempty(trialTypeFilter)
    %work out logic of trial type selection
    logStr = '>=|<=|==|~=|>|<'; %these are the characters allowed
    logIdx_start = regexp(trialTypeFilter,logStr,'start');
    logIdx_end = regexp(trialTypeFilter,logStr,'end');
    SelField = trialTypeFilter(1:logIdx_start-1);
    SelLogic = trialTypeFilter(logIdx_start:logIdx_end);
    SelVal = trialTypeFilter((logIdx_end+1):end);
    %make the selection
    if strcmp(SelField(end-4:end),'times')
        delay = trialsT.(SelField) - trialsT.stimOn_times; %for now, assume we want the time relative to stimOn
        trialSel = eval(['delay',SelLogic,SelVal])';
    else
        trialSel = eval(['trialsT.(SelField)',SelLogic,SelVal])';
    end
else
    trialSel = true(1,size(trialsT,1));
end

%if no trial type values were defined, take all unique values
if isempty(trialTypeVals) && isfield(trialsT,trialTypeField)
    trialTypeVals = unique(trialsT.(trialTypeField));
end

%now make our final selection
trialCnd = logical([]);
if length(trialTypeVals)>1
    for iV = 1:length(trialTypeVals)
        trialCnd(iV,:) = trialSel & trialsT.(trialTypeField)'==trialTypeVals(iV);
    end
else
    trialCnd = trialSel;
end
