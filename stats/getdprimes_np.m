function dprimes = getdprimes_np(allResps, trialCnd)

% computes 'non-parametric dprimes' of each neuron's avg response across two sets of
% trials, using the mann-whitney U and concerting this to a pseudo d-prime
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

if isa(trialCnd,'double')
    trialCnd = logical(trialCnd);
end

nNeurons = size(allResps,1);
nSignal = sum(trialCnd(2,:));
nNoise = sum(trialCnd(1,:));

% Stack signal and noise horizontally
combined = [allResps(:,trialCnd(2,:)), allResps(:,trialCnd(1,:))]; % [nComparisons x (nSignal+nNoise)]

% Compute ranks row-wise
ranks = zeros(size(combined));
for i = 1:nNeurons
    ranks(i,:) = tiedrank(combined(i,:));
end

% Sum of ranks for signal samples
rankSumSignal = sum(ranks(:,1:nSignal),2); % sum across signal columns

% Compute Mann-Whitney U
U = rankSumSignal - nSignal*(nSignal+1)/2;

% Compute AUC
AUC = U ./ (nSignal * nNoise);

% clip AUC away from 0 and 1 with a scale-aware epsilon
eps = 1/(2*nSignal*nNoise); 
AUC = min(max(AUC, eps), 1-eps);

% Convert AUC to d-prime
dprimes = sqrt(2) * erfinv(2*AUC - 1);
