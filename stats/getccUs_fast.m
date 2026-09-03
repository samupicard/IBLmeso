function Us = getccUs_fast(resps, trialCondition, trialVariable)
% resps: nROIs x nTrials
% trialCondition: 1 x nTrials (or nTrials x 1) integer labels
% trialVariable:  1 x nTrials logical (A vs B)

trialCondition = trialCondition(:)'; 
trialVariable  = trialVariable(:)';

nROIs = size(resps,1);
uCond = unique(trialCondition);

% Precompute indices once (only ~6 conditions)
condIdx  = cell(numel(uCond),1);
condPosA = cell(numel(uCond),1);
nA = zeros(numel(uCond),1);
nB = zeros(numel(uCond),1);
validCond = false(numel(uCond),1);

for ci = 1:numel(uCond)
    idx = find(trialCondition == uCond(ci));
    condIdx{ci} = idx;

    posA = find(trialVariable(idx));
    condPosA{ci} = posA;

    nA(ci) = numel(posA);
    nB(ci) = numel(idx) - nA(ci);
    validCond(ci) = (nA(ci) > 0) && (nB(ci) > 0);
end

dTotal = sum(nA(validCond) .* nB(validCond));
if dTotal == 0
    Us = nan(1,nROIs);
    return
end

nTotal = zeros(nROIs,1);

% Chunk ROIs to reduce tiedrank cost
block = 512; % tune: 256/512/1024 depending on CPU/cache

for r0 = 1:block:nROIs
    r1 = min(nROIs, r0+block-1);

    nBlock = zeros(r1-r0+1,1);

    for ci = 1:numel(uCond)
        if ~validCond(ci), continue; end

        idx = condIdx{ci};
        posA = condPosA{ci};  % positions within idx

        X = single(resps(r0:r1, idx));   % blockROIs x nTcond

        % tiedrank ranks along columns, so transpose to nTcond x blockROIs
        R = tiedrank(X');      % nTcond x blockROIs
        % sum ranks of A trials
        R1 = sum(R(posA, :), 1)';  % blockROIs x 1

        numer = R1 - nA(ci) * (nA(ci)+1) / 2;
        nBlock = nBlock + numer;
    end

    nTotal(r0:r1) = nBlock;
end

Us = (nTotal ./ dTotal).';
end
