
function [psychoM,RTs] = getPsycho_fitBlock(trialsT,blocks,ks)
%GETPSYCHO_FITBLOCK Psychometric fit with delayed subjective block transitions.
%
% [psychoM,RTs] = getPsycho_fitBlock(trialsT,blocks)
%x
% blocks : true (default) to fit block-dependent choice bias.
%
% The full psychometric model is:
% P(right) = logistic(v*sign(X)*abs(X)^q + db*B + b)
%
% B is the subjective block prior, shifted by k trials relative
% to the objective probabilityLeft sequence.
%
% Positive k = delayed block switch.
% Negative k = early block switch.
%
% The best lag maximizes the incremental deviance explained of the full
% model, relative to a null model without bias blocks.
%
% NB: Descriptive psychometric and RT statistics are subsequently
% grouped by subjective rather than objective block probability.

if nargin < 2
    blocks = true;
    ks = 4; %this is the average best lag across most expert mice
end
if nargin < 3
    ks = 4; %this is the average best lag across most expert mice
end

%% Load trial variables

% Classify engaged and NoGo trials
trialsT = classifyEngagementTrials(trialsT);

% Trials eligible for psychometric analysis
validEngagedGo = ~trialsT.disengagedTrials & ~trialsT.noGoTrials;

% get block probabilities
trialsT.probabilityLeft = round(trialsT.probabilityLeft,1);
probaLeft = trialsT.probabilityLeft(:);
if ~blocks
    probaLeft(:) = 0.5;
end
probas = unique(probaLeft(~isnan(probaLeft)));
if numel(probas) > 3 % Handle trainingCW debiasing probabilities
    probaLeft(:) = 0.5;
    probas = 0.5;
end
[~,i] = sort(abs(probas-0.5),'ascend'); % Place 0.5 first, if present
probas = probas(i);

choice = -trialsT.choice(:);
choiceTypes = unique(choice(~isnan(choice) & choice~=0));
numChoiceTypes = numel(choiceTypes);

rt = trialsT.response_times(:) - trialsT.stimOn_times(:);

if ismember('contrastDiff',trialsT.Properties.VariableNames)
    cDiff = trialsT.contrastDiff(:);
else
    cRight = trialsT.contrastRight(:);
    cLeft = trialsT.contrastLeft(:);
    cRight(isnan(cRight)) = 0;
    cLeft(isnan(cLeft)) = 0;
    cDiff = cRight-cLeft;
end

%% Select contrasts with sufficient repeats

cVals = unique(cDiff(~isnan(cDiff)));

cKeep = true(size(cVals));
for i = 1:ceil(numel(cVals)/2)
    if sum(cDiff==cVals(i))<4 || sum(cDiff==-cVals(i))<4
        cKeep(i) = false;
        cKeep(end+1-i) = false;
    end
end
cVals = cVals(cKeep);

%% Fit psychometric model and estimate subjective prior lag

F = 'logit';

psychoM.bestLag = NaN;
psychoM.lagScan = [];
psychoM.subjectiveProbLeft = probaLeft;
psychoM.devExplained = NaN;
psychoM.devExplainedBlock = NaN;
psychoM.fracDevExplainedBlock = NaN;

if numel(probas)>1

    % Objective block regressor:
    % low pLeft = +1 (rightward bias)
    % high pLeft = -1 (leftward bias)
    B = zeros(size(probaLeft));
    B(probaLeft==probas(end-1)) =  1;
    B(probaLeft==probas(end))   = -1;

    dbs = nan(size(ks));
    devs = nan(size(ks));
    fits = cell(size(ks));
    stats = cell(size(ks));

    % Build lagged regressors independently within each session
    Bscan = zeros(numel(B),numel(ks));

    for ik = 1:numel(ks)
        idx = getSessionLagIndices(trialsT,ks(ik));
        Bscan(:,ik) = B(idx);
    end

    % Same fitting trials for every lag, making deviances comparable
    validTrials = validEngagedGo & ...
        round(probaLeft,1)~=0.5 & ...
        all(Bscan~=0,2) & ...
        ismember(choice,choiceTypes) & ...
        ~isnan(cDiff);

    assert(numChoiceTypes==2, ...
        'Expected two choice types for psychometric fitting.');

    y = double(choice(validTrials)==choiceTypes(2));
    X = cDiff(validTrials);

    for ik = 1:numel(ks)
        Bk = Bscan(validTrials,ik);

        [fits{ik},devs(ik),stats{ik}] = ...
            fitPM_bias(y,X,Bk,F);

        dbs(ik) = fits{ik}.db;
    end

    % Null model (intercept only)
    [~,nullDev] = glmfit(ones(size(y)),y, ...
        'binomial','link',F,'constant','off');

    % Reduced model (contrast only)
    [~,devNoBlock] = fitPM_nobias(y,X,F);

    % Deviance explained
    devExplained = 1-devs/nullDev;
    devExplainedBlock = (devNoBlock-devs)/nullDev;

    % Select lag with strongest incremental DE and expected bias direction
    validLags = dbs > 0 & isfinite(devExplainedBlock);

    if any(validLags)
        score = devExplainedBlock;
        score(~validLags) = -Inf;
        [~,iBest] = max(score);
        bestLag = ks(iBest);
    else
        if length(ks)>1
            iBest = find(ks==3); %default to k=3
        else
            iBest = 1;
        end
        bestLag = ks(iBest);
    end

    
    %Refit selected model using all eligible trials for this lag
    Bk = Bscan(:,iBest);
    validTrials = validEngagedGo & ...
        round(probaLeft,1)~=0.5 & ...
        Bk~=0 & ...
        ismember(choice,choiceTypes) & ...
        ~isnan(cDiff);
    y = double(choice(validTrials)==choiceTypes(2));
    X = cDiff(validTrials);
    [coeff,dev,stat] = fitPM_bias(y,X,Bk(validTrials),F);

    %Recompute deviance explained for the final model
    [~,nullDevFinal] = glmfit(ones(size(y)),y, 'binomial','Link',F,'constant','off');
    [~,devNoBlockFinal] = fitPM_nobias(y,X,F);
    
    % Store different DE metrics for the final model
    psychoM.devExplained = 1 - dev/nullDevFinal;
    psychoM.devExplainedBlock = (devNoBlockFinal-dev)/nullDevFinal; 
    psychoM.fracDevExplainedBlock = (devNoBlockFinal - dev) / (nullDevFinal - dev);

    % Store lag scan table
    psychoM.bestLag = bestLag;
    psychoM.lagScan = table(ks(:),dbs(:),devs(:), ...
        devExplained(:),devExplainedBlock(:), ...
        'VariableNames',{'k','deltabias','deviance', ...
        'devExplained','devExplainedBlock'});

    % Shift block probabilities independently within each session
    idx = getSessionLagIndices(trialsT,bestLag);
    psychoM.subjectiveProbLeft = probaLeft(idx);

else

    % No block-dependent bias
    validTrials = validEngagedGo & ismember(choice,choiceTypes) & ~isnan(cDiff);
    y = double(choice(validTrials)==choiceTypes(end));
    X = cDiff(validTrials);

    [coeff,dev,stat] = fitPM_nobias(y,X,F);

end

%% Descriptive psychometric and RT statistics
% Group trials by SUBJECTIVE block probabilities

subjProbaLeft = psychoM.subjectiveProbLeft;

nBlocks = numel(probas);
nContrasts = numel(cVals);

psychoM.blockProbLeft = probas;
psychoM.contrastDiffs = cVals;
psychoM.choiceTypes = choiceTypes;

psychoM.numTrials = zeros(nBlocks,nContrasts);
psychoM.numChooseR = zeros(nBlocks,numChoiceTypes,nContrasts);
psychoM.mean = nan(nBlocks,numChoiceTypes,nContrasts);
psychoM.ci = nan(nBlocks,numChoiceTypes,nContrasts);

RTs.blockProbLeft = probas;
RTs.contrastDiffs = cVals;
RTs.choiceTypes = choiceTypes;
RTs.mean = nan(nBlocks,numChoiceTypes,nContrasts);
RTs.CIlo = nan(nBlocks,numChoiceTypes,nContrasts);
RTs.CIhi = nan(nBlocks,numChoiceTypes,nContrasts);

for iBlock = 1:nBlocks
    for iContrast = 1:nContrasts

        incl = validTrials & ...
            cDiff==cVals(iContrast) & ...
            subjProbaLeft==probas(iBlock);

        n = sum(incl);
        psychoM.numTrials(iBlock,iContrast) = n;

        if n==0
            continue
        end

        for iChoice = 1:numChoiceTypes

            nChoice = sum(incl & choice==choiceTypes(iChoice));
            p = nChoice/n;

            psychoM.numChooseR(iBlock,iChoice,iContrast) = nChoice;
            psychoM.mean(iBlock,iChoice,iContrast) = p;
            psychoM.ci(iBlock,iChoice,iContrast) = ...
                1.96*sqrt(p*(1-p)/n);

            % RT statistics for this choice and contrast
            thisRT = rt(incl & choice==choiceTypes(iChoice));
            thisRT = thisRT(~isnan(thisRT));

            if ~isempty(thisRT)
                qRT = quantile(thisRT,[0.25 0.5 0.75]);
                RTs.mean(iBlock,iChoice,iContrast) = qRT(2);
                RTs.CIlo(iBlock,iChoice,iContrast) = qRT(2)-qRT(1);
                RTs.CIhi(iBlock,iChoice,iContrast) = qRT(3)-qRT(2);
            end
        end
    end
end

%% Store fitted parameters

psychoM.bias = coeff.b;
psychoM.sensitivity = coeff.v;
psychoM.exponent = coeff.q;

if numel(probas)>1
    psychoM.deltabias = coeff.db;
else
    psychoM.deltabias = NaN;
end

psychoM.dev = dev;
psychoM.stat = stat;

end


function idx = getSessionLagIndices(trialsT,k)
% Shift trial indices by k, clipping at each session boundary.

n = height(trialsT);
idx = (1:n)';

if ismember('sessionID',trialsT.Properties.VariableNames)
    sessionID = trialsT.sessionID;
else
    sessionID = ones(n,1);
end

for s = unique(sessionID(:))'
    ii = find(sessionID==s);
    idx(ii) = ii(max(1,min(numel(ii),(1:numel(ii))'-k)));
end

end