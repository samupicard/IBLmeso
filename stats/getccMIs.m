function MIs = getccMIs(resps, trialsT, trialVariable, trialTypeField, condFields, minTrialN)
% Vectorized Modulation Index: computes condition-combined Modulation Index for all ROIs
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

miTot = zeros(nROIs,1);
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

    mi = (mRespsA - mRespsB) ./ (mRespsA + mRespsB);

    miTot = miTot + mi;
    cnt = cnt+1;
end

if cnt == 0
    MIs = nan(1,nROIs);
else
    MIs = miTot / cnt;
end

end