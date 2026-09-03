function meandiffs = getccMeanDiffs(resps, trialsT, trialVariable, trialTypeField, condFields, minTrialN)
% Vectorized Mean Diff: computes condition-combined Mean Differences for all ROIs
% resps: nROIs x nTrials
% trialVariable: 1 x nTrials logical (or nTrials x 1)
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

trialVariable = trialVariable(:)';  % force row logical

% Condition labels per trial
trialsT_conds = trialsT(:, conds_to_match);
[~, ~, trialCondition] = unique(trialsT_conds, 'rows');
trialCondition = trialCondition(:)'; % row

uCond = unique(trialCondition);
nROIs = size(resps,1);

diffTot = zeros(nROIs,1);
cnt = 0;

for ci = 1:numel(uCond)
    inclT = (trialCondition == uCond(ci));

    A = trialVariable & inclT;
    B = (~trialVariable) & inclT;

    nA = sum(A);
    nB = sum(B);
    if nA<minTrialN || nB<minTrialN
        continue
    end

    mRespsA = mean(resps(:,A), 2, 'omitnan');
    mRespsB = mean(resps(:,B), 2, 'omitnan');

    diffTot = diffTot + (mRespsA - mRespsB);
    cnt = cnt+1;
end

if cnt == 0
    meandiffs = nan(1,nROIs);
else
    meandiffs = diffTot / cnt;
end

end