function dprimes = getdprimes(allResps, trialCnd, K, minTrialN)

% computes dprimes of each neuron's avg response across two sets of trials.
% optionally uses k-fold cross-validation
%
%inputs
%   allResps [neurons x trials] avg response of each neuron to each trial
%   trialCnd [2 x trials] logical array of trials to compare w/ each other
%   NB: computes resps in LAST set in trialCnd minus resps in FIRST set in trialCnd.

%outputs
%   dprimes [1, neurons]
%
%
% written by Samuel Picard (Oct 2023)
% update Dec 2025: added cross-validation option

if isa(trialCnd,'double')
    trialCnd = logical(trialCnd);
end

% Indices for the two conditions:
% First row = "noise" (condition 1), last row = "signal" (condition 2)
noiseIdx  = trialCnd(1,:);      % FIRST set
signalIdx = trialCnd(end,:);    % LAST set

% Only consider trials that belong to either of the two sets
usedTrials = noiseIdx | signalIdx;

noiseTrials  = find(noiseIdx);
signalTrials = find(signalIdx);

nNoise  = numel(noiseTrials);
nSignal = numel(signalTrials);

if nargin<4
    minTrialN = 5;  % Default minimum number of trials if not provided
end

% Basic sanity checks
if nNoise < minTrialN || nSignal < minTrialN
    warning('Not enough trials in one or both conditions. Returning NaNs.');
    dprimes = nan(1, size(allResps,1));
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

    % Randomly permute trials within each condition and split into K folds
    noiseTrials  = noiseTrials(randperm(nNoise));
    signalTrials = signalTrials(randperm(nSignal));
    
    % Fold boundaries
    noiseEdges  = round(linspace(0, nNoise,  K+1));
    signalEdges = round(linspace(0, nSignal, K+1));

    d_fold = nan(nNeurons, K);

    % ==============================
    % PER-FOLD CROSS-VALIDATION LOOP
    % ==============================
    for k = 1:K

        % ---- Select test trials ----
        testNoise  = noiseTrials( noiseEdges(k)+1  : noiseEdges(k+1) );
        testSignal = signalTrials(signalEdges(k)+1 : signalEdges(k+1) );

        testIdx = false(1, size(allResps,2));
        testIdx(testNoise)  = true;
        testIdx(testSignal) = true;

        % ---- Training trials ----
        trainNoiseIdx  = noiseIdx  & ~testIdx;
        trainSignalIdx = signalIdx & ~testIdx;

        % Extract training responses
        respNoiseTrain  = allResps(:, trainNoiseIdx);   % [neurons x nTrainNoise]
        respSignalTrain = allResps(:, trainSignalIdx);  % [neurons x nTrainSignal]

        % Compute μ and σ on training data (vectorized)
        mu1 = nanmean(respNoiseTrain,  2);   % noise mean
        mu2 = nanmean(respSignalTrain, 2);   % signal mean

        sigma1 = nanstd(respNoiseTrain,  0, 2);   % noise SD
        sigma2 = nanstd(respSignalTrain, 0, 2);   % signal SD

        % Pooled SD
        sp = sqrt(0.5 * (sigma1.^2 + sigma2.^2));

        % Avoid divide-by-zero
        invalid = (sp == 0) | isnan(sp);
        sp(invalid) = NaN;

        % Classical d' formula for this fold
        d_fold(:, k) = (mu2 - mu1) ./ sp;
    end

    % Average across folds
    dprimes_col = nanmean(d_fold, 2);

    % Return as 1 x neurons
    dprimes = dprimes_col.';

else
    difs = nanmean(allResps(:,signalIdx),2) - nanmean(allResps(:,noiseIdx ),2);
    vars = sqrt(0.5*(var(allResps(:,noiseIdx),0,2,'omitnan') + var(allResps(:,signalIdx),0,2,'omitnan')));
    dprimes = difs ./ vars;
end

% Neurons with no usable trials end up with d' = NaN
dprimes(~isfinite(dprimes)) = NaN;