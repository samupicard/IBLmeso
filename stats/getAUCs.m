function AUCs = getAUCs(allResps, trialCnd, K)

% computes AUCs of each neuron's avg response across two sets of trials.
% optionally uses k-fold cross-validation
%
%inputs
%   allResps [neurons x trials] avg response of each neuron to each trial
%   trialCnd [2 x trials] logical array of trials to compare w/ each other
%
%outputs
%   AUCs [1, neurons]
%
%NB: AUC is the Mann-Whitney U statistic scaled by the total number of
%condition-pairs.
%
% written by Samuel Picard (Dec 2025)

if isa(trialCnd,'double')
    trialCnd = logical(trialCnd);
end

% Indices for the two conditions:
% First row = "noise" (condition 1), last row = "signal" (condition 2)
noiseIdx  = trialCnd(1,:);      % FIRST set
signalIdx = trialCnd(end,:);    % LAST set

noiseTrials  = find(noiseIdx);
signalTrials = find(signalIdx);

if any(noiseIdx & signalIdx)
    error('trialCnd rows must be mutually exclusive.');
end

nNoise  = numel(noiseTrials);
nSignal = numel(signalTrials);

% Basic sanity checks
if nNoise < 2 || nSignal < 2
    warning('Not enough trials in one or both conditions. Returning NaNs.');
    AUCs = nan(1, size(allResps,1));
    return;
end

if nargin<3 || isempty(K)
    K = false;
elseif islogical(K) && K
    K = min(10, min(nNoise, nSignal));
elseif K > 1
    K = min(K, min(nNoise, nSignal));
elseif ~K
    K = false;
else
    warning('argument K is invalid. Returning NaNs.');
    dprimes = nan(1, size(allResps,1));
    return;
end

if K
    
    if K < 2
        warning('K < 2 after adjustment. Returning NaNs.');
        dprimes = nan(1, size(allResps,1));
        return;
    end
    
    [nNeurons, nTrials] = size(allResps);

    
    % Randomize trials within each condition
    noiseTrials  = noiseTrials(randperm(nNoise));
    signalTrials = signalTrials(randperm(nSignal));

    % CV fold boundaries
    noiseEdges  = round(linspace(0, nNoise,  K+1));
    signalEdges = round(linspace(0, nSignal, K+1));

    % Store AUC for each fold
    AUC_folds = nan(nNeurons, K);

    % ==============================
    % Cross-validation loop
    % ==============================
    for k = 1:K

        % Test trials for this fold
        testNoise  = noiseTrials( noiseEdges(k)+1  : noiseEdges(k+1) );
        testSignal = signalTrials(signalEdges(k)+1 : signalEdges(k+1) );

        testIdx = false(1, nTrials);
        testIdx(testNoise)  = true;
        testIdx(testSignal) = true;

        respTest = allResps(:, testIdx);     % [neurons x nTestTrials]

        % Labels for test trials: 0=noise, 1=signal
        labels = zeros(1, sum(testIdx));
        labels(ismember(find(testIdx), testSignal)) = 1;

        % Separate responses within fold
        noiseResp  = respTest(:, labels == 0);
        signalResp = respTest(:, labels == 1);

        nNoiseT  = size(noiseResp, 2);
        nSignalT = size(signalResp, 2);

        if nNoiseT == 0 || nSignalT == 0
            continue;
        end

        % Combine and rank
        combined = [noiseResp signalResp];           % [neurons x total]
        ranks    = tiedrank(combined.').';           % rank each neuron's row

        % Sum of ranks of signal responses
        R_signal = sum(ranks(:, nNoiseT+1:end), 2);

        % Mann-Whitney U
        U = R_signal - nSignalT*(nSignalT+1)/2;

        % AUC per neuron
        AUC = U ./ (nSignalT * nNoiseT);

        AUC_folds(:, k) = AUC;
    end

    % Average across folds, return as [1 x neurons]
    AUCs = nanmean(AUC_folds, 2).';

else
        noiseResp  = allResps(:, noiseIdx);    % [neurons × nNoise]
    signalResp = allResps(:, signalIdx);   % [neurons × nSignal]

    nNoise  = size(noiseResp, 2);
    nSignal = size(signalResp, 2);

    if nNoise == 0 || nSignal == 0
        AUCs = nan(1, size(allResps,1));
        return;
    end

    % Combine responses for ranking
    combined = [noiseResp signalResp];   % [neurons x totalTrials]

    % Rank responses across all trials per neuron
    ranks = tiedrank(combined.').';      % rank rows using transpose trick

    % Sum of ranks of signal responses
    R_signal = sum(ranks(:, nNoise+1:end), 2);

    % Mann–Whitney U statistic
    U = R_signal - nSignal*(nSignal + 1)/2;

    % Convert to AUC
    AUCs = (U ./ (nSignal * nNoise)).';

    % OPTIONAL: Convert to unsigned AUC
    %AUCs = max(AUCs, 1 - AUCs);
    
end
    