function [psychoM,RTs] = getPsycho_fit(trialsT,blocks)

%this function gets the psychometric fit data from a trials table.
%by default, separates by block type (if there are blocks)
%run getPsycho_fit(trialsT, false) to NOT separate by block type.
%
% if there are blocks, uses glm model fitting to get bias difference between blocks

if nargin<2
    blocks = true;
end

if blocks
    probaLeft = trialsT.probabilityLeft;
    probas = unique(probaLeft);
else %pretend that everything is equiprobable
    probaLeft = 0.5*ones(size(trialsT.probabilityLeft));
    probas = 0.5;
end

%quick patch to deal with debiasing probabilities in trainingCW
if length(probas)>3
    probaLeft = 0.5*ones(size(probaLeft));
    probas = 0.5;
end

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

% Fit PM data using this function: y = 1/(1+exp(-(v*X^q+db*B+b))); 
% which is a logistic function with delta-bias 'db' and overall bias 'b'

% Fit psychometric model
F = 'logit';

psychoM.bestLag = NaN;
psychoM.lagScan = [];

if length(probas)>1

    % Objective block regressor across ALL trials
    % (including 50% trials, to preserve trial indexing)
    B = zeros(size(probaLeft));
    B(probaLeft==probas(end-1)) = 1;
    B(probaLeft==probas(end))   = -1;

    % Candidate subjective block-switch lags (in trials)
    ks = -5:15;

    dbs  = nan(size(ks));
    devs = nan(size(ks));
    fits = cell(size(ks));
    stats = cell(size(ks));

    for ik = 1:numel(ks)
        k = ks(ik);

        % Shift B by k trials, padding with boundary values
        idx = max(1, min(numel(B), (1:numel(B)) - k));
        Bk = B(idx);

        % Exclude objective 50% trials and undefined subjective priors
        validTrials = round(probaLeft,1)~=0.5 & Bk~=0 ...
            & ~isnan(choice) & ~isnan(cDiff);

        y = choice(validTrials)==choiceTypes(2);
        X = cDiff(validTrials);

        [fits{ik},devs(ik),stats{ik}] = ...
            fitPM_bias(y,X,Bk(validTrials),F);

        dbs(ik) = fits{ik}.db;
    end

    % Fit reduced model (contrast only, no block prior)
    [~,devNoBlock] = fitPM_nobias(y,X,F);

    % Compute deviance explained relative to intercept-only model
    validTrials = round(probaLeft,1)~=0.5 & ...
        ~isnan(choice) & ~isnan(cDiff);
    y = choice(validTrials)==choiceTypes(2);
    [~,nullDev] = glmfit(ones(sum(validTrials),1),y, ...
        'binomial','constant','off');

    % Deviance explained by block regressor, relative to a null model
    devExplained = 1 - devs/nullDev;

    % Incremental deviance explained by block regressor
    devExplainedBlock = (devNoBlock - devs) / nullDev;

    % Select best lag 
    %[~,iBest] = max(dbs); %maximizing delta-bias
    [~,iBest] = max(devExplainedBlock); %maximizing incremental deviance explained


    coeff = fits{iBest};
    dev   = devs(iBest);
    stat  = stats{iBest};

    psychoM.bestLag = ks(iBest);


    psychoM.lagScan = table(ks(:),dbs(:),devs(:), ...
        devExplained(:),devExplainedBlock(:), ...
        'VariableNames',{'k','deltabias','deviance', ...
        'devExplained','devExplainedBlock'});

else

    y = choice==choiceTypes(2);
    X = cDiff;

    [coeff,dev,stat] = fitPM_nobias(y,X,F);

end

psychoM.bias = coeff.b;
if length(probas)>1
    psychoM.deltabias = coeff.db;
else
    psychoM.deltabias = NaN;
end
psychoM.sensitivity = coeff.v;
psychoM.exponent = coeff.q;
psychoM.dev = dev;
psychoM.stat = stat;


end

