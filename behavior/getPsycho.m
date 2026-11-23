function [psychoM,RTs] = getPsycho(trialsT)
%this function gets the psychometric data from a trials table (ignoring block)

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

%get all the relevant data
psychoM = {};
RTs = {};
psychoM.contrastDiffs = cVals;
psychoM.choiceTypes = choiceTypes;
psychoM.numTrials = zeros(1,length(cVals));
psychoM.numChooseR = zeros(numChoiceTypes, length(cVals));
psychoM.mean = zeros(numChoiceTypes,length(cVals));
psychoM.CI = zeros(numChoiceTypes,length(cVals));

RTs.contrastDiffs = cVals;
RTs.choiceTypes = choiceTypes;
RTs.mean = zeros(numChoiceTypes,length(cVals));
RTs.CIlo = zeros(numChoiceTypes,length(cVals));
RTs.CIhi = zeros(numChoiceTypes,length(cVals));

for iChoice = 1:numChoiceTypes
    for iContrast = 1:length(cVals)
        incl = cDiff==cVals(iContrast); %Logical array of trials of this contrast + blockType
        psychoM.numTrials(iContrast) = max([1,sum(incl)]);
        psychoM.numChooseR(iChoice,iContrast) = sum(choice==choiceTypes(iChoice)&incl);
        psychoM.mean(iChoice, iContrast) = psychoM.numChooseR(iChoice,iContrast)/psychoM.numTrials(iContrast);
        psychoM.ci(iChoice, iContrast) = 1.96*sqrt(psychoM.mean(iChoice, iContrast)*(1-psychoM.mean(iChoice, iContrast))/psychoM.numTrials(iContrast));
        
        q = quantile(rt(choice==choiceTypes(iChoice)&incl), 3);
        RTs.mean(iChoice, iContrast) = q(2);
        RTs.CIlo(iChoice, iContrast) = q(2)-q(1);
        RTs.CIhi(iChoice, iContrast) = q(3)-q(2);
    end
end


%calculate the params for psycho.meanetric line fit
xx = cVals'*100; %Unique contrasts in percent
nn = psychoM.numTrials;
if size(psychoM.mean,1)>1
    pp = squeeze(psychoM.mean(2,:)); %Performance
else
    pp = squeeze(psychoM.mean(:,:));
end
if any(xx==0)
    intercept = psychoM.mean(xx==0);
else
    intercept=0;
end

%Calculate params for line fit using Matteo's function
%syntax: mle_fit_psycho(data,P_model,parstart,parmin,parmax,nfits), where
%order of params = bias, 1/slope (threshold yielding 80%), gamma1, gamma2
psychoM.pars = mle_fit_psycho([xx;nn;pp],'erf_psycho_2gammas',...
    [mean(xx) 20 min(pp) 1-max(pp)],...
    [min(xx) 0.5*min(diff(xx)), 0, 0],...
    [max(xx), 50, .5, .5],25); %same but with two different lapse rates and more sensible start, min and max pars.
