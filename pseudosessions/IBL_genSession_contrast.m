function pseudoSession = IBL_genSession_contrast(sz, probabilityLeft)

%use this to generate new trial contrasts according to ephysChoiceWorld heuristics
%randomizing the contrasts, but NOT the blocks

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

pseudoSession = {};
pseudoSession.probabilityLeft = probabilityLeft(1:sz);
pseudoSession.contrastLeft = contrastLeft;
pseudoSession.contrastRight = contrastRight;

end