
function results = comparePsychometricModelsCV(trialsT, varargin)
%COMPAREPSYCHOMETRICMODELSCV Blocked cross-validation of psychometric priors.
%
% results = comparePsychometricModelsCV(trialsT, ...
%     'lags',1:20,'taus',1:2:100,'nFolds',5,'qs',0.05:0.05:1);
%
% Models:
%   baseline : logistic(b + v * sign(X)*abs(X)^q)
%   block    : baseline + db * lagged block (+1 for pLeft=.2)
%   side     : baseline + db * exponentially averaged previous sides
%   contrast : baseline + db * exponentially averaged transformed contrasts
%   drift    : logistic(b0 + v*f_q(X) + latent bias deviation)
%              sigma estimated on training choices by Laplace evidence;
%              online predictions assimilate each test choice ONLY AFTER
%              predicting it (optional strict held-out mode).
%              Requires psychometricDriftCV.m
%
% q is selected ONCE from the contrast-only baseline on each outer training
% fold, and shared by all models in that fold. A full-data q is also stored.
% If 'fixedQ' is supplied, the externally chosen value is used throughout.
% History calculations reset at original session boundaries and use all
% observed stimuli, including those on trials excluded from the GLM.
% Tau is in nonzero-stimulus presentations (zero contrasts do not update
% or decay history). History at trial t excludes stimulus t.
%
% Identical eligibility for all models and folds: engaged Go, valid choice,
% finite contrast, objective pLeft=.2/.8, and nonzero lagged B at ALL
% candidate lags. Thus biased trials shortly after initial .5 periods
% may be removed, regardless of which model is evaluated.
%
% This function assumes classifyEngagementTrials.m is on the MATLAB path.
%
% Outputs:
%   results.cvTable          one row per model, held-out deviance and deltaDE
%   results.folds            held-out predictions, foldwise q, k, taus, sigma
%   results.full             best full-data fits and tuning choices
%   results.validTrials      logical mask into original trialsT
%   results.foldId           fold for each original trial (0 = excluded)

ip = inputParser;
ip.addParameter('lags',1:20,@(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x==round(x)));
ip.addParameter('taus',1:2:100,@(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)&x>0));
ip.addParameter('nFolds',5,@(x)isscalar(x)&&x>=2&&x==round(x));
ip.addParameter('qs',0.05:0.05:1,@(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)&x>0));
ip.addParameter('fixedQ',[],@(x)isempty(x)||(isscalar(x)&&isfinite(x)&&x>0));
ip.addParameter('sigmaBounds',[.002 .5],@(x)isnumeric(x)&&numel(x)==2&&all(x>0)&&x(1)<x(2));
ip.addParameter('initialBiasSD',.5,@(x)isscalar(x)&&isfinite(x)&&x>0);
ip.addParameter('nQuadrature',21,@(x)isscalar(x)&&x>=5&&x==round(x));
ip.addParameter('driftEvaluation','online', ...
    @(x)any(strcmpi(char(x),{'online','strict'})));
ip.parse(varargin{:});
opt = ip.Results;
ks = unique(opt.lags(:)','stable');
taus = unique(opt.taus(:)','stable');
qs = unique(opt.qs(:)','stable');
nFolds = opt.nFolds;

% Recompute eligibility rather than trusting potentially stale columns.
trialsT = classifyEngagementTrials(trialsT);
n = height(trialsT);
if ismember('contrastDiff',trialsT.Properties.VariableNames)
    X = double(trialsT.contrastDiff(:));
else
    cr = double(trialsT.contrastRight(:));
    cl = double(trialsT.contrastLeft(:));
    cr(isnan(cr))=0;
    cl(isnan(cl))=0;
    X = cr-cl;
end
choice = -double(trialsT.choice(:)); % +1 is right after IBL conversion
y = double(choice==1);
pLeft = double(trialsT.probabilityLeft(:));

% A session ID change means a reset of stimulus history and shifted blocks.
if ismember('sessionID',trialsT.Properties.VariableNames)
    sessionCode = findgroups(trialsT.sessionID);
elseif ismember('sessionNr',trialsT.Properties.VariableNames)
    sessionCode = findgroups(trialsT.sessionNr);
else
    sessionCode = ones(n,1);
end
sessionCode = sessionCode(:);
assert(all(~isnan(sessionCode)), 'Missing session identifiers.');
starts = [1; find(diff(sessionCode)~=0)+1];
ends = [starts(2:end)-1; n];

% Objective block sign: +1 indicates right-biased (pLeft=.2).
B = zeros(n,1);
B(abs(pLeft-.2)<1e-6) = 1;
B(abs(pLeft-.8)<1e-6) = -1;
Bscan = zeros(n,numel(ks));
for i = 1:numel(starts)
    ii = (starts(i):ends(i))';
    for j = 1:numel(ks)
        localIdx = max(1,min(numel(ii),(1:numel(ii))'-ks(j)));
        Bscan(ii,j) = B(ii(localIdx));
    end
end

% Include all engaged Go trials with valid choices and contrasts,
% regardless of objective or subjective block probability
valid = ~trialsT.noGoTrials(:) & ...
        ~trialsT.disengagedTrials(:) & ...
        ismember(choice,[-1 1]) & ...
        isfinite(X);
validIdx = find(valid);
nValid = numel(validIdx);
assert(nValid>=nFolds*10,'Too few eligible trials for requested CV folds.');

% Assign contiguous blocks in the sequence of eligible trials.
foldId = zeros(n,1);
cut = floor((0:nFolds)*nValid/nFolds);
for f = 1:nFolds
    foldId(validIdx(cut(f)+1:cut(f+1))) = f;
end

names = {'baseline','block','side','contrast','drift'};
nModels = numel(names);
cvPrediction = nan(n,nModels);
foldDeviance = nan(nFolds,nModels);
foldQ = nan(nFolds,1);
foldBestK = nan(nFolds,1);
foldBestTauSide = nan(nFolds,1);
foldBestTauContrast = nan(nFolds,1);
foldDb = nan(nFolds,nModels);
foldSigma = nan(nFolds,1);
foldSigmaAtBoundary = false(nFolds,1);
foldBeta = cell(nFolds,nModels);
cvBiasPredict = nan(n,1);
trainCriterion = nan(nFolds,nModels);

for f = 1:nFolds
    test = valid & foldId==f;
    train = valid & foldId~=f;
    assert(numel(unique(y(train)))==2,'Training fold has only one choice class.');

    % Fit q once using baseline training data; never test choices.
    if isempty(opt.fixedQ)
        q = selectQ(X(train),y(train),qs);
    else
        q = opt.fixedQ;
    end
    foldQ(f) = q;
    fX = signedPower(X,q);

    % Common contrast-only baseline.
    betaBase = glmfit(fX(train),y(train),'binomial','link','logit');
    foldBeta{f,1} = betaBase;
    cvPrediction(test,1) = logistic([ones(nnz(test),1), fX(test)]*betaBase);
    trainCriterion(f,1) = bernoulliDev(y(train), ...
        logistic([ones(nnz(train),1), fX(train)]*betaBase));

    % Generate history regressors using all stimuli, not merely fit trials.
    Hside = nan(n,numel(taus));
    Hcontrast = nan(n,numel(taus));
    for j = 1:numel(taus)
        [Hside(:,j),Hcontrast(:,j)] = ...
            historyPriors(X,taus(j),q,starts,ends);
    end

    matrices = {Bscan,Hside,Hcontrast};
    for m = 2:4
        H = matrices{m-1};
        devTrain = nan(1,size(H,2));
        betas = cell(1,size(H,2));
        for j = 1:size(H,2)
            design = [fX(train),H(train,j)];
            [betas{j},devTrain(j)] = ...
                glmfit(design,y(train),'binomial','link','logit');
        end
        [trainCriterion(f,m),jBest] = min(devTrain);
        beta = betas{jBest};
        foldBeta{f,m} = beta;
        foldDb(f,m) = beta(3);
        if m==2, foldBestK(f) = ks(jBest); end
        if m==3, foldBestTauSide(f) = taus(jBest); end
        if m==4, foldBestTauContrast(f) = taus(jBest); end
        cvPrediction(test,m) = logistic( ...
            [ones(nnz(test),1),fX(test),H(test,jBest)]*beta);
    end

    % Fit sigma using ONLY training choices. In online mode, each held-out
    % choice updates the filter AFTER its own prediction. Predictions
    % contain no information from the current or later choice.
    if strcmpi(opt.driftEvaluation,'online')
        driftUpdate = valid;
    else
        driftUpdate = train;
    end
    drift = psychometricDriftCV(betaBase(1)+betaBase(2)*fX, ...
        y,train,starts,ends, ...
        'sigmaBounds',opt.sigmaBounds, ...
        'initialBiasSD',opt.initialBiasSD, ...
        'nQuadrature',opt.nQuadrature, ...
        'updateMask',driftUpdate);
    foldBeta{f,5} = betaBase; 
    cvBiasPredict(test) = drift.biasPredict(test); % One-step-ahead latent bias deviation for held-out trials
    foldSigma(f) = drift.sigma;
    foldSigmaAtBoundary(f) = drift.sigmaAtBoundary;
    cvPrediction(test,5) = drift.pPredict(test);
    trainCriterion(f,5) = -2*drift.logEvidence; % Not a training deviance.

    for m = 1:nModels
        foldDeviance(f,m) = bernoulliDev(y(test),cvPrediction(test,m));
    end
end

% An intercept-only predictive null, fitted on each training fold.
nullPrediction = nan(n,1);
for f = 1:nFolds
    test = valid & foldId==f;
    train = valid & foldId~=f;
    pTrain = mean(y(train));
    nullPrediction(test) = pTrain;
end
nullDev = bernoulliDev(y(valid),nullPrediction(valid));
cvDev = sum(foldDeviance,1);
cvDeltaDE = (cvDev(1)-cvDev)/nullDev;
cvDevExplained = 1-cvDev/nullDev;

results.cvTable = table(string(names(:)),cvDev(:),cvDevExplained(:), ...
    cvDeltaDE(:),'VariableNames', ...
    {'model','cvDeviance','cvDevExplained','cvDeltaDE'});
results.validTrials = valid;
results.foldId = foldId;
results.folds = table((1:nFolds)',foldQ,foldBestK,foldBestTauSide, ...
    foldBestTauContrast,foldSigma,foldSigmaAtBoundary,'VariableNames', ...
    {'fold','q','bestK','bestTauSide','bestTauContrast', ...
     'sigmaDrift','sigmaDriftAtBoundary'});
results.folds.deviance = foldDeviance;
results.folds.db = foldDb;
results.folds.trainDeviance = trainCriterion;
% For drift, the fifth column of trainDeviance is -2 log evidence,
% not deviance. See results.folds.driftNeg2LogEvidence.
results.folds.trainDeviance(:,5) = NaN;
results.folds.driftNeg2LogEvidence = trainCriterion(:,5);
results.cvPrediction = cvPrediction;
results.cvNullPrediction = nullPrediction;
results.cvNullDeviance = nullDev;
results.foldDeviance = foldDeviance;
results.driftEvaluation = char(opt.driftEvaluation);
results.lags = ks;
results.taus = taus;
results.foldBeta = foldBeta;
results.cvBiasPredict = cvBiasPredict;

% Fit a final descriptive model using all eligible trials, with one q.
if isempty(opt.fixedQ)
    qFull = selectQ(X(valid),y(valid),qs);
else
    qFull = opt.fixedQ;
end
fX = signedPower(X,qFull);
results.full.q = qFull;
[results.full.baseline.beta,results.full.baseline.deviance] = ...
    glmfit(fX(valid),y(valid),'binomial','link','logit');
[Hside,Hcontrast] = precomputeHistories(X,taus,qFull,starts,ends);
matrices = {Bscan,Hside,Hcontrast};
for m=2:4
    H = matrices{m-1};
    devs = nan(1,size(H,2));
    betas = cell(1,size(H,2));
    for j=1:size(H,2)
        [betas{j},devs(j)] = glmfit([fX(valid),H(valid,j)], ...
            y(valid),'binomial','link','logit');
    end
    [bestDev,jBest] = min(devs);
    model = names{m};
    results.full.(model).beta = betas{jBest};
    results.full.(model).deviance = bestDev;
    results.full.(model).db = betas{jBest}(3);
    results.full.(model).candidateDeviance = devs;
    if m==2
        results.full.(model).bestK = ks(jBest);
    else
        results.full.(model).bestTau = taus(jBest);
    end
end
driftFull = psychometricDriftCV( ...
    results.full.baseline.beta(1)+ ...
    results.full.baseline.beta(2)*fX,y,valid,starts,ends, ...
    'sigmaBounds',opt.sigmaBounds, ...
    'initialBiasSD',opt.initialBiasSD, ...
    'nQuadrature',opt.nQuadrature);
results.full.drift.sigma = driftFull.sigma;
results.full.drift.logEvidence = driftFull.logEvidence;
results.full.drift.sigmaAtBoundary = driftFull.sigmaAtBoundary;
results.full.drift.pPredict = driftFull.pPredict;
results.full.drift.biasPredict = driftFull.biasPredict;
results.full.drift.biasPredictSD = driftFull.biasPredictSD;
results.full.drift.biasFilter = driftFull.biasFilter;
results.full.drift.biasFilterSD = driftFull.biasFilterSD;
results.full.drift.initialBiasSD = opt.initialBiasSD;
results.full.drift.v = results.full.baseline.beta(2);
results.full.drift.b0 = results.full.baseline.beta(1);
results.full.nValidTrials = nValid;
end

function q = selectQ(X,y,qs)
% Choose q using only contrast and an intercept.
deviance = nan(numel(qs),1);
for j=1:numel(qs)
    [~,deviance(j)] = glmfit(signedPower(X,qs(j)),y, ...
        'binomial','link','logit');
end
[~,idx] = min(deviance);
q = qs(idx);
end

function z = signedPower(x,q)
z = sign(x).*abs(x).^q;
end

function p = logistic(z)
p = 1./(1+exp(-z));
p = max(1e-12,min(1-1e-12,p));
end

function d = bernoulliDev(y,p)
p = max(1e-12,min(1-1e-12,p));
d = -2*sum(y.*log(p)+(1-y).*log(1-p));
end

function [side,contrast] = historyPriors(X,tau,q,starts,ends)
% Decay is per nonzero contrast presentation; output excludes current X.
alpha = exp(-1/tau);
side = zeros(size(X));
contrast = zeros(size(X));
for s = 1:numel(starts)
    numS=0; denS=0; numC=0; denC=0;
    for t = starts(s):ends(s)
        if denS>0, side(t)=numS/denS; end
        if denC>0, contrast(t)=numC/denC; end
        x = X(t);
        if isfinite(x) && x~=0
            w = abs(x)^q;
            numS = alpha*numS + sign(x);
            denS = alpha*denS + 1;
            numC = alpha*numC + sign(x)*w;
            denC = alpha*denC + w;
        end
    end
end
end

function [Hside,Hcontrast] = precomputeHistories(X,taus,q,starts,ends)
Hside = zeros(numel(X),numel(taus));
Hcontrast = zeros(numel(X),numel(taus));
for j=1:numel(taus)
    [Hside(:,j),Hcontrast(:,j)] = ...
        historyPriors(X,taus(j),q,starts,ends);
end
end
