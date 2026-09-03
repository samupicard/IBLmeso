function pseudoSession = IBL_genSession(sz, trialsT, block_len_lambda, block_len_min, block_len_max)

%generate new blocks and trials according to ephysChoiceWorld heuristics
%update Dec 2023: this is now compatible with multiple trial bouts
%update May 2024: include full trials Table

if nargin < 3
    block_len_lambda=60;
    block_len_min=20;
    block_len_max=100;
    %initblock_len=90; %initial 50/50 block (this is now inferred form pLeft)
end

block_len = uint16(round(rnd.exp(block_len_lambda,[1 round(size(trialsT,1)/block_len_min)],[block_len_min,block_len_max])));

pLeft = trialsT.probabilityLeft;

if size(pLeft,1)>1
    pLeft=pLeft';
end

%find initial 50/50 blocks
initblock_edges = [1,diff(pLeft==0.5)];
initblocks(1,:) = find(initblock_edges==1); 
initblocks(2,:) = find(initblock_edges==-1);
boutsizes = diff([initblocks(1,:) sz+1]);
initblock_lens = diff(initblocks);

%generate block probabilities (keeping 50/50 blocks the same)
probabilityLeft = [0.5*ones(1,sz)];
iBlock = 1;
for iBout = 1:length(boutsizes)
    pLeft_bout = [0.5*ones(1,boutsizes(iBout))];
    iTrial = initblock_lens(iBout)+1;
    block_probs = [0.2 0.8];
    curr_block = randsample(block_probs,1);
    while iTrial<=boutsizes(iBout)
        if curr_block == block_probs(1)
            curr_block = block_probs(2);
        elseif curr_block == block_probs(2)
            curr_block = block_probs(1);
        end
        pLeft_bout(iTrial:iTrial+block_len(iBlock))=curr_block;
        iTrial = iTrial+block_len(iBlock)+1;
        iBlock = iBlock+1;
    end
    probabilityLeft(initblocks(1,iBout)+(0:boutsizes(iBout)-1)) = pLeft_bout(1:boutsizes(iBout));
end

%then use this to generate new trial contrasts
contrastVals = [0 0.0625 0.125 0.25 1];
trialSides = ['L','R'];
contrastLeft = nan(1,sz);
contrastRight = nan(1,sz);
for iTrial = 1:sz
    side = randsample(trialSides,1,true,[probabilityLeft(iTrial) 1-probabilityLeft(iTrial)]);
    if strcmp(side,'L')
        contrastLeft(iTrial) = randsample(contrastVals,1);
    elseif strcmp(side,'R')
        contrastRight(iTrial) = randsample(contrastVals,1);
    end
end
contrastDiff = getContrastDiff(contrastLeft, contrastRight);

%add prior
if 1
    tau = 10;
    prior = [];
    if any(ismember(trialsT.Properties.VariableNames,'sessionNr'))
        sessionNr = trialsT{:,'sessionNr'};
    else
        sessionNr = ones(1,size(trialsT,1));
    end
    sessionNr_unique = unique(sessionNr);
    for iS = 1:length(sessionNr_unique)
        iTrials = sessionNr==sessionNr_unique(iS);
        isLeftTrial = -sign(contrastDiff(iTrials));
        isLeftTrial(contrastDiff(iTrials)==0) = nan;
        isLeftTrial(isLeftTrial==-1) = 0;
        prior = [prior; computeExponentialBias(isLeftTrial,tau)];
    end
end

pseudoSession = {};
pseudoSession.probabilityLeft = probabilityLeft';
%pseudoSession.contrastLeft = contrastLeft;
%pseudoSession.contrastRight = contrastRight;
pseudoSession.contrastDiff = [getContrastDiff(contrastLeft,contrastRight)]';
pseudoSession.choice = trialsT.choice;
pseudoSession.feedbackType = trialsT.feedbackType;
%pseudoSession.prior = prior;

% pseudoSession = table( ...
%     probabilityLeft', ...
%     contrastDiff', ...
%     trialsT.choice, ...
%     trialsT.feedbackType, ...
%     prior, ...
%     'VariableNames', {'probabilityLeft', 'contrastDiff', 'choice', 'feedbackType', 'prior'});

end