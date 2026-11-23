function results = getPsycho_fit(trialsT,taus)
%GETPSYCHO_FIT Fit PM models with stimulus-history priors.
%
% results = getPsycho_fit(trialsT,taus)
%
% Compares two continuous stimulus-history regressors:
%   'side'     : exponentially averaged stimulus side
%   'contrast' : stimulus side weighted by absolute contrast
%
% Selects tau by maximizing incremental deviance explained.
%
% Reuses fitPM_bias and fitPM_nobias.

if nargin<2 || isempty(taus)
    taus = [1:2:39 45:5:100];
end

taus = taus(:);
assert(all(isfinite(taus) & taus>0), ...
    'All decay constants must be positive and finite.');

%% Load trial variables

trialsT = classifyEngagementTrials(trialsT);

choice = -trialsT.choice(:);

if ismember('contrastDiff',trialsT.Properties.VariableNames)
    cDiff = trialsT.contrastDiff(:);
else
    cRight = trialsT.contrastRight(:);
    cLeft = trialsT.contrastLeft(:);
    cRight(isnan(cRight)) = 0;
    cLeft(isnan(cLeft)) = 0;
    cDiff = cRight-cLeft;
end

% Session identifiers
if ismember('sessionID',trialsT.Properties.VariableNames)
    sessionID = trialsT.sessionID;
else
    sessionID = ones(height(trialsT),1);
end

%% Select valid trials

% Include 0.5 blocks: history is independent of block identity.
% Do not remove excluded trials from the stimulus history.

validTrials = ~trialsT.noGoTrials & ...
              ~trialsT.disengagedTrials & ...
              ismember(choice,[-1 1]) & ...
              isfinite(cDiff);

y = double(choice(validTrials)==1);
X = cDiff(validTrials);

assert(numel(unique(y))==2, ...
    'Valid trials must contain both left and right choices.');

%% Fit null and contrast-only models

F = 'logit';

[~,nullDev] = glmfit(ones(size(y)),y, ...
    'binomial','link',F,'constant','off');

[coeffNoHistory,devNoHistory,statNoHistory] = ...
    fitPM_nobias(y,X,F);

%% Fit history models

modes = {'side','contrast'};

results = struct();
results.taus = taus;
results.validTrials = validTrials;
results.nullDev = nullDev;
results.devNoHistory = devNoHistory;
results.coeffNoHistory = coeffNoHistory;
results.statNoHistory = statNoHistory;

for iMode = 1:numel(modes)

    mode = modes{iMode};

    nTau = numel(taus);

    dbs = nan(nTau,1);
    devs = nan(nTau,1);
    fits = cell(nTau,1);
    stats = cell(nTau,1);

    for iTau = 1:nTau

        tau = taus(iTau);

        % Compute stimulus history over all trials
        H = computeStimulusHistory( ...
            cDiff,tau,sessionID,mode);

        % Select same trials as reduced/null models
        Hfit = H(validTrials);

        % Fit model using existing GLM helper
        [fits{iTau},devs(iTau),stats{iTau}] = ...
            fitPM_bias(y,X,Hfit,F);

        dbs(iTau) = fits{iTau}.db;
    end

    %% Incremental deviance explained

    devExplained = 1-devs/nullDev;

    devExplainedHistory = ...
        (devNoHistory-devs)/nullDev;

    %% Select optimal decay constant

    [~,iBest] = max(devExplainedHistory);

    bestTau = taus(iBest);

    % Reconstruct best history regressor
    bestHistory = computeStimulusHistory( ...
        cDiff,bestTau,sessionID,mode);

    %% Store results

    M = struct();

    M.mode = mode;
    M.bestTau = bestTau;
    M.history = bestHistory;

    M.bias = fits{iBest}.b;
    M.sensitivity = fits{iBest}.v;
    M.exponent = fits{iBest}.q;
    M.deltabias = fits{iBest}.db;

    M.dev = devs(iBest);
    M.stat = stats{iBest};

    M.devExplained = devExplained(iBest);
    M.devExplainedHistory = devExplainedHistory(iBest);

    if nullDev>devs(iBest)
        M.fracDevExplainedHistory = ...
            (devNoHistory-devs(iBest))/ ...
            (nullDev-devs(iBest));
    else
        M.fracDevExplainedHistory = NaN;
    end

    M.tauScan = table(taus,dbs,devs, ...
        devExplained,devExplainedHistory, ...
        'VariableNames',{'tau','deltabias','deviance', ...
        'devExplained','devExplainedHistory'});

    results.(mode) = M;

end

end
