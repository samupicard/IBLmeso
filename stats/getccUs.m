function Us = getccUs(resps, trialsT, trialVariable, trialTypeField, condFields, minTrialN)
% Vectorized ccU: computes condition-combined Mann-Whitney U for all ROIs
% resps: nROIs x nTrials
% trialVariable: 1 x nTrials logical (or nTrials x 1)
% trialsT: table with trial metadata
%
% Samuel Picard / optimized suggestion

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

nTotal = zeros(nROIs,1);
dTotal = 0;

for ci = 1:numel(uCond)
    inclT = (trialCondition == uCond(ci));

    A = trialVariable & inclT;
    B = (~trialVariable) & inclT;

    nA = sum(A);
    nB = sum(B);
    if nA<minTrialN || nB<minTrialN
        continue
    end

    idx = find(inclT);        % trials in this condition
    %idxA = find(A(inclT));    % positions within idx (not global)
    relA = trialVariable(inclT);

    % responses for this condition: nROIs x nT
    X = resps(:, idx);

    % rank trials within each ROI (column-wise in tiedrank)
    R = tiedrank(X');   % nT x nROIs
    R = R';             % nROIs x nT

    % sum ranks for group A (all ROIs at once)
    %R1 = sum(R(:, idxA), 2);
    R1 = sum(R(:, relA), 2);

    % Mann-Whitney U numerator for each ROI
    numer = R1 - nA*(nA+1)/2;
    %numer = (nA*nB) - (R1 - nA*(nA+1)/2); %flipped test direction

    %actual_AUC_this_loop = mean(numer ./ (nA*nB)); 
    %printf('Cond %d: nA=%d, nB=%d, Mean AUC=%.3f\n', ci, nA, nB, actual_AUC_this_loop);

    nTotal = nTotal + numer;
    dTotal = dTotal + nA*nB;
end

if dTotal == 0
    Us = nan(1,nROIs);
else
    Us = (nTotal ./ dTotal).';  % row: 1 x nROIs
end
end