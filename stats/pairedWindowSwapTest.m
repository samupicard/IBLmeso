function [observedStat, nullStats, pValues, significant, nValidTrials] = ...
    pairedWindowSwapTest( ...
        pairedDifference, ...
        nPermutations, ...
        nTrialsMin, ...
        pthresh)
%
% pairedWindowSwapTest
%
% Randomly swapping response and baseline within a trial changes:
%
%   response - baseline
%
% to:
%
%   baseline - response
%
% Therefore each permutation is implemented as an independent random sign
% flip of each trial's paired difference.

nROIs = size(pairedDifference, 1);
nTrials = size(pairedDifference, 2);

validValues = isfinite(pairedDifference);

nValidTrials = sum(validValues, 2);

observedStat = mean( ...
    pairedDifference, ...
    2, ...
    'omitnan');

observedStat(nValidTrials < nTrialsMin) = nan;

nullStats = nan(nROIs, nPermutations);

% Generate the random trial swaps once. A value of -1 means response and
% baseline are swapped for that trial.
swapSigns = 2 * (rand(nTrials, nPermutations) >= 0.5) - 1;

for iPermutation = 1:nPermutations

    signedDifference = ...
        pairedDifference .* swapSigns(:, iPermutation)';

    nullStats(:, iPermutation) = mean( ...
        signedDifference, ...
        2, ...
        'omitnan');
end

nullStats(nValidTrials < nTrialsMin, :) = nan;

pValues = nan(nROIs, 1);

for iROI = 1:nROIs

    if ~isfinite(observedStat(iROI))
        continue
    end

    nullThisROI = nullStats(iROI, :);
    nullThisROI = nullThisROI(isfinite(nullThisROI));

    if isempty(nullThisROI)
        continue
    end

    % Two-sided randomization p-value with the standard +1 correction.
    pValues(iROI) = ...
        (1 + sum( ...
            abs(nullThisROI) >= abs(observedStat(iROI)))) / ...
        (numel(nullThisROI) + 1);
end

significant = pValues < pthresh;

end

