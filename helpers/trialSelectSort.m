function [trialSortIdx, trialCnd_sorted] = trialSelectSort(trialsT, trialTypeField, trialTypeVals, trialTypeFilter, evnt_lock, evnt_sort)

%selects and sorts trials based on a number of conditions
%
% Samuel Picard

%append event name with _times
s = '_times'; N = min(length(evnt_lock),length(s));
if strcmp(evnt_lock(end-(N-1):end), s)
    evnt_nm = evnt_lock;
    evnt_lock = evnt_nm(1:end-length(s));
else
    evnt_nm = [evnt_lock s];
end

if strcmp(evnt_sort, 'type')
sortType = true;
end    


evnt_sort = 'type'; %if set to 'type', will group by trialTypeVals. If set to some timed event (e.g. feedback) will sort by latency relative to time-locking event. If empty, no sorting.
sortType = true;

%make trial selection by condition
trialCnd = trialSelection(trialsT,trialTypeField,trialTypeVals,trialTypeFilter);

%get latencies to sort trials by
sortLats = [];
if (ischar(evnt_sort) || isstring(evnt_sort)) && ~isempty(evnt_sort) && ~strcmp(evnt_sort,'type')
    if any(contains(trialsT.Properties.VariableNames,[evnt_sort '_times']))
        sortLats = trialsT.([evnt_sort,'_times']) -  trialsT.(evnt_nm);
    end
end


if sortType && ~isempty(evnt_sort) && any(contains(trialsT.Properties.VariableNames,[evnt_sort '_times']))
    trialSortIdx = [];
    [~,subSortIdx] = sort(abs(sortLats),'ascend');
    trialCnd_presorted = trialCnd(:,subSortIdx);
    for iCnd = 1:length(trialTypeVals)
        trials_of_this_type = find(trialCnd(iCnd,:));
        [~,isort] = sort(abs(sortLats(trials_of_this_type)),'ascend');
        trialSortIdx = [trialSortIdx trials_of_this_type(isort)];
    end
    trialCnd_sorted = trialCnd(:,trialSortIdx);
else
    if strcmp(evnt_sort,'type')
        trialSortIdx = [];
        for iCnd = 1:length(trialTypeVals)
            trialSortIdx = [trialSortIdx find(trialCnd(iCnd,:))];
        end
    elseif ~isempty(evnt_sort) && any(contains(trialsT.Properties.VariableNames,[evnt_sort '_times']))
        [~,trialSortIdx] = sort(sortLats,'ascend');
    else
        trialSortIdx = 1:size(trialCnd,2);
    end
    trialCnd_sorted = trialCnd(:,trialSortIdx);
end
    