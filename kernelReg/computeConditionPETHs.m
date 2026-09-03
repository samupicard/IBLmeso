function peth = computeConditionPETHs(P, Y, predictorInfo, options)
%COMPUTECONDITIONPETHS Compute condition-specific event-triggered averages.
%
% Y can be either:
%   - observed neural activity
%   - model-predicted neural activity
%
% P is time × predictors.
% Y is time × neurons.
%
%   predictorInfo
%       Table containing:
%           column
%           name
%           type
%           lagSeconds
%
% Output
%   peth
%       Struct array with one element per discrete-lagged predictor:
%
%       peth(ii).name
%       peth(ii).lags
%       peth(ii).mean
%       peth(ii).sem
%       peth(ii).nEvents
%       peth(ii).conditionValues
%       peth(ii).columns
%       peth(ii).zeroLagColumn
%
%       mean and sem are neurons × lags.
%
% Options
%   eventThreshold
%       A design-matrix value is considered active when its absolute
%       value exceeds this threshold.
%
%   baselineWindow
%       Optional two-element window, in seconds, used to baseline-subtract
%       each neuron's PETH. For example:
%
%           baselineWindow = [-1 -0.5]
%
%       Default is empty, meaning no baseline subtraction.
%
% Notes
%   This assumes that each discrete-lagged column is a shifted version of
%   an event indicator, with zero values outside event occurrences.

arguments
    P double
    Y double
    predictorInfo table

    options.eventThreshold (1,1) double = 0
    options.baselineWindow double = double.empty
end

%% Basic checks

if size(P, 1) ~= size(Y, 1)
    error( ...
        "P and Y must have the same number of rows. " + ...
        "P has %d rows and Y has %d rows.", ...
        size(P, 1), size(Y, 1));
end

if max(predictorInfo.column) > size(P, 2)
    error( ...
        "predictorInfo references column %d, but P has only %d columns.", ...
        max(predictorInfo.column), size(P, 2));
end

%% Find discrete-lagged predictors

isLagged = ...
    predictorInfo.type == "discrete-lagged" & ...
    isfinite(predictorInfo.lagSeconds);

laggedInfo = predictorInfo(isLagged, :);

predictorNames = unique(laggedInfo.name, "stable");
nPredictors = numel(predictorNames);
nNeurons = size(Y, 2);

if nPredictors == 0
    error("No discrete-lagged predictors were found.");
end

%% Allocate output

peth = struct( ...
    "name", cell(nPredictors, 1), ...
    "lags", cell(nPredictors, 1), ...
    "mean", cell(nPredictors, 1), ...
    "sem", cell(nPredictors, 1), ...
    "nEvents", cell(nPredictors, 1), ...
    "conditionValues", cell(nPredictors, 1), ...
    "columns", cell(nPredictors, 1), ...
    "zeroLagColumn", cell(nPredictors, 1));

%% Compute PETH for each predictor

for ii = 1:nPredictors
    predictorName = predictorNames(ii);

    groupMask = laggedInfo.name == predictorName;

    columns = laggedInfo.column(groupMask);
    lags = laggedInfo.lagSeconds(groupMask);

    % Put lag bins in chronological order.
    [lags, lagOrder] = sort(lags);
    columns = columns(lagOrder);

    nLags = numel(lags);

    zeroLagColumn = columns( ...
    find(abs(lags) == min(abs(lags)), 1));

eventValues = P(:, zeroLagColumn);

conditionValues = unique( ...
    eventValues(isfinite(eventValues) & eventValues ~= 0));

conditionValues = conditionValues(:).';
nConditions = numel(conditionValues);

meanPETH = nan(nNeurons, nLags, nConditions);
semPETH = nan(nNeurons, nLags, nConditions);
nEvents = zeros(nLags, nConditions);

for conditionIndex = 1:nConditions
    conditionValue = conditionValues(conditionIndex);

    for lagIndex = 1:nLags
        predictorColumn = P(:, columns(lagIndex));

        % Allow a small tolerance for floating-point contrast values.
        tolerance = max(1e-10, abs(conditionValue) * 1e-8);

        eventMask = ...
            isfinite(predictorColumn) & ...
            abs(predictorColumn - conditionValue) <= tolerance;

        eventActivity = Y(eventMask, :);

        nEvents(lagIndex, conditionIndex) = sum(eventMask);

        if isempty(eventActivity)
            continue
        end

        meanPETH(:, lagIndex, conditionIndex) = ...
            mean(eventActivity, 1, "omitnan").';

        nFinite = sum(isfinite(eventActivity), 1);
        activitySD = std(eventActivity, 0, 1, "omitnan");

        thisSEM = activitySD ./ sqrt(nFinite);
        thisSEM(nFinite == 0) = NaN;

        semPETH(:, lagIndex, conditionIndex) = thisSEM.';
    end
end

    %% Optional baseline subtraction

    if ~isempty(options.baselineWindow)
        if numel(options.baselineWindow) ~= 2
            error("baselineWindow must contain [startTime endTime].");
        end

        baselineMask = ...
            lags >= options.baselineWindow(1) & ...
            lags <= options.baselineWindow(2);

        if ~any(baselineMask)
            error( ...
                "No lag bins fall inside baseline window [%.3g %.3g] " + ...
                "for predictor '%s'.", ...
                options.baselineWindow(1), ...
                options.baselineWindow(2), ...
                predictorName);
        end

        baseline = mean(meanPETH(:, baselineMask), 2, "omitnan");
        meanPETH = meanPETH - baseline;
    end

    %% Identify the column closest to zero lag

    [~, zeroLagIndex] = min(abs(lags));
    zeroLagColumn = columns(zeroLagIndex);

    %% Store results

    peth(ii).name = predictorName;
    peth(ii).lags = lags(:).';
    peth(ii).mean = meanPETH;
    peth(ii).sem = semPETH;
    peth(ii).nEvents = nEvents;
    peth(ii).conditionValues = conditionValues;
    peth(ii).columns = columns(:).';
    peth(ii).zeroLagColumn = zeroLagColumn;
end

end