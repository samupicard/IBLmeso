function [psychoM,RTs] = getPsycho_erffit(trialsT)

%this function gets the psychometric fit data from a trials table (separating
%by block type)
%
% use glm model fitting to get bias difference between blocks

probaLeft = trialsT.probabilityLeft;
probas = unique(probaLeft);

[~,i] = sort(abs(probas-0.49),'ascend');
probas = probas(i); %quick fix to make sure we have .5 first (if it's there)

choice = -trialsT.choice;
feedback = trialsT.feedbackType;
choiceTypes = unique(choice(choice~=0 * isnan(choice))); % -1 is RIGHT (IBL)
numChoiceTypes = numel(choiceTypes);

rt = trialsT.response_times - trialsT.stimOn_times;

if any(strcmp('contrastDiff',trialsT.Properties.VariableNames))
    cDiff = trialsT.contrastDiff;
else
    trials_contrastRight = trialsT.contrastRight;
    trials_contrastRight(isnan(trials_contrastRight))=0;
    trials_contrastLeft = trialsT.contrastLeft;
    trials_contrastLeft(isnan(trials_contrastLeft))=0;
    cDiff = trials_contrastRight-trials_contrastLeft;
end
cVals = unique(cDiff);

%take out those contrasts for which there are fewer than 4 repeats
cKeep = true(1,length(cVals));
for i=1:ceil(length(cVals)/2)
    if (sum(cDiff==cVals(i))<4) || (sum(cDiff==-cVals(i))<4)
        cKeep(i)=false;
        cKeep(end+1-i)=false;
    end
end
cVals=cVals(cKeep);

%take out 0.5 constrasts (only appear in training)
%cVals = cVals(abs(cVals)~=0.5);

%get general descriptive stats (divided by block)
psychoM = {};
RTs = {};
psychoM.blockProbLeft = probas;
psychoM.contrastDiffs = cVals;
psychoM.choiceTypes = choiceTypes;
psychoM.numTrials = zeros(length(probas),length(cVals));
psychoM.numChooseR = zeros(length(probas),numChoiceTypes, length(cVals));
psychoM.mean = zeros(length(probas),numChoiceTypes,length(cVals));
psychoM.CI = zeros(length(probas),numChoiceTypes,length(cVals));

RTs.blockProbLeft = probas;
RTs.contrastDiffs = cVals;
RTs.choiceTypes = choiceTypes;
RTs.mean = zeros(length(probas),numChoiceTypes,length(cVals));
RTs.CIlo = zeros(length(probas),numChoiceTypes,length(cVals));
RTs.CIhi = zeros(length(probas),numChoiceTypes,length(cVals));

for iBlock = 1:length(probas)
    for iChoice = 1:numChoiceTypes
        for iContrast = 1:length(cVals)
            incl = cDiff==cVals(iContrast)&probaLeft==probas(iBlock); %Logical array of trials of this contrast + blockType
            psychoM.numTrials(iBlock,iContrast) = max([1,sum(incl)]);
            psychoM.numChooseR(iBlock,iChoice,iContrast) = sum(choice==choiceTypes(iChoice)&incl);
            psychoM.mean(iBlock, iChoice, iContrast) = psychoM.numChooseR(iBlock,iChoice,iContrast)/psychoM.numTrials(iBlock,iContrast);
            psychoM.ci(iBlock, iChoice, iContrast) = 1.96*sqrt(psychoM.mean(iBlock,iChoice, iContrast)*(1-psychoM.mean(iBlock,iChoice, iContrast))/psychoM.numTrials(iBlock,iContrast));
            
            q = quantile(rt(choice==choiceTypes(iChoice)&incl), 3);
            RTs.mean(iBlock, iChoice, iContrast) = q(2);
            RTs.CIlo(iBlock, iChoice, iContrast) = q(2)-q(1);
            RTs.CIhi(iBlock, iChoice, iContrast) = q(3)-q(2);
        end
    end
end

% Fit PM data using old function

xx = cVals';
all_pars = [];
for iP=[1,3]
    nn = psychoM.numTrials(iP,:);
    pp = squeeze(psychoM.mean(iP,1,:))';
    pars = mle_fit_psycho([xx;nn;pp],'erf_psycho_2gammas',[mean(xx) 20 min(pp) 1-max(pp)],[min(xx) 0.5*min(diff(xx)), 0, 0],[max(xx), 50, .5, .5],25); %same but with two different lapse rates and more sensible start, min and max pars.
    if iP==1, all_pars = pars; all_nn = nn; else all_pars = [all_pars; pars]; all_nn = [all_nn; nn]; end
end

psychoM.pars = all_pars;

end

