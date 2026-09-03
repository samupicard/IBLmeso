function means = getccMeans(resps, trialsT, trialTypeField, condFields, minTrialN)
% Vectorized Mean Diff: computes condition-combined Means for all ROIs
% resps: nROIs x nTrials
% trialsT: table with trial metadata
%
% Samuel Picard

if nargin<6
    minTrialN = 5;  % Default minimum number of trials if not provided
end
if nargin<5
    condFields = '';
end

if isempty(condFields)
    condFields = getBalancedConds(trialTypeField);
end
conds_to_match = condFields(~strcmp(condFields, trialTypeField));

% Condition labels per trial
trialsT_conds = trialsT(:, conds_to_match);
[~, ~, trialCondition] = unique(trialsT_conds, 'rows');
trialCondition = trialCondition(:)'; % row

uCond = unique(trialCondition);
nROIs = size(resps,1);

Tot = zeros(nROIs,1);
cnt = 0;

for ci = 1:numel(uCond)
    inclT = (trialCondition == uCond(ci));

    n = sum(inclT);
    if n<minTrialN
        continue
    end

    mResps = mean(resps(:,inclT), 2, 'omitnan');

    Tot = Tot + mResps;
    cnt = cnt+1;
end

if cnt == 0
    means = nan(1,nROIs);
else
    means = Tot / cnt;
end

end