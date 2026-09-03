function plotFullReducedModelPETHs( ...
    P, F, predictorInfo, results, neuronIndices, options)
%PLOTFULLREDUCEDMODELPETHS Compare full- and reduced-model PETH predictions.
%
% One figure is produced per neuron. Columns correspond to predictors.
%
%   First row:
%       Empirical mean, solid line
%       Empirical SEM, shaded region
%
%   Second row:
%       Full-model prediction, solid line
%
%   Third row:
%       Reduced-model prediction, dashed line
%
% For each predictor p, the reduced model is:
%
%       results.predictor(p).otherModel
%
% which was fitted using all design-matrix columns except those belonging
% to predictor p.
%
% INPUTS
% ------
% P
%   Full T-by-D design matrix used to fit the full model.
%
% F
%   T-by-N measured neural activity.
%
% predictorInfo
%   Metadata describing the columns of P. This is passed to
%   computeConditionPETHs.
%
% results
%   Output from fitIBLKernelSelectivity. Expected fields include:
%
%       results.full.model
%       results.predictor(p).name
%       results.predictor(p).otherModel
%       results.predictor(p).otherColumns
%
% neuronIndices
%   Indices of neurons to plot.
%
% OPTIONS
% -------
% predictorNames
%   Predictors to display. By default, all predictors with an available
%   reduced model are used.
%
% showSEM
%   Display empirical SEM shading. Default: true.
%
% showEmpiricalMean
%   Display the empirical mean as a faint dotted line. Default: false.
%
% showLegend
%   Display condition legends. Default: false.
%
% baselineSubtract
%   Baseline-subtract empirical and predicted PETHs. Default: false.
%
% baselineWindow
%   Baseline interval in seconds. Default: [-1 -0.5].
%
% fullLineWidth
%   Width of full-model prediction. Default: 1.75.
%
% reducedLineWidth
%   Width of reduced-model prediction. Default: 1.75.
%
% empiricalLineWidth
%   Width of optional empirical mean. Default: 1.
%
% semAlpha
%   Transparency of empirical SEM bounds. Default: 0.15.
%
% shareYAcrossRows
%   Use the same y-limits in both rows for each neuron. Default: true.

arguments
    P double
    F double
    predictorInfo
    results struct
    neuronIndices (:,1) double {mustBeInteger,mustBePositive}

    options.predictorNames string = strings(0,1)

    options.showSEM (1,1) logical = true
    options.showEmpiricalMean (1,1) logical = false
    options.showLegend (1,1) logical = false

    options.baselineSubtract (1,1) logical = false
    options.baselineWindow (1,2) double = [-1 -0.5]

    options.fullLineWidth (1,1) double {mustBePositive} = 1.75
    options.reducedLineWidth (1,1) double {mustBePositive} = 1.75
    options.empiricalLineWidth (1,1) double {mustBePositive} = 1

    options.semAlpha (1,1) double ...
        {mustBeGreaterThanOrEqual(options.semAlpha,0), ...
        mustBeLessThanOrEqual(options.semAlpha,1)} = 0.2

    options.shareYAcrossRows (1,1) logical = true
end

%% Validate dimensions

nFrames = size(P,1);
nNeurons = size(F,2);

if size(F,1) ~= nFrames
    error( ...
        "plotFullReducedModelPETHs:FrameCountMismatch", ...
        "P and F must have the same number of rows.");
end

if any(neuronIndices > nNeurons)
    error( ...
        "plotFullReducedModelPETHs:NeuronIndex", ...
        "A requested neuron index exceeds the available %d neurons.", ...
        nNeurons);
end

if ~isfield(results,"full") || ...
        ~isfield(results.full,"model") || ...
        isempty(results.full.model)

    error( ...
        "plotFullReducedModelPETHs:MissingFullModel", ...
        "results.full.model is missing or empty.");
end

if ~isfield(results,"predictor") || isempty(results.predictor)
    error( ...
        "plotFullReducedModelPETHs:MissingReducedModels", ...
        "results.predictor is missing or empty.");
end

%% Determine available discrete-lagged predictors

resultPredictorNames = string({results.predictor.name});

hasReducedModel = arrayfun( ...
    @(x) isfield(x,"otherModel") && ~isempty(x.otherModel), ...
    results.predictor);

% Identify predictors whose predictorInfo type is "discrete-lagged".
predictorInfoNames = string(results.predictorInfo.name);
predictorInfoTypes = string(results.predictorInfo.type);

isDiscreteLagged = false(size(resultPredictorNames));

for predictorIndex = 1:numel(resultPredictorNames)

    matchingRows = ...
        predictorInfoNames == resultPredictorNames(predictorIndex);

    if any(matchingRows)
        isDiscreteLagged(predictorIndex) = all( ...
            strcmpi( ...
            predictorInfoTypes(matchingRows), ...
            "discrete-lagged"));
    end
end

availablePredictor = hasReducedModel' & isDiscreteLagged;
availablePredictorNames = ...
    resultPredictorNames(availablePredictor);

if isempty(options.predictorNames)

    % By default, plot all discrete-lagged predictors that have a fitted
    % reduced model.
    predictorNames = availablePredictorNames;

else

    % Explicitly requested predictors must also be discrete-lagged and
    % have an available reduced model.
    requestedPredictorNames = options.predictorNames(:).';

    invalidPredictorNames = requestedPredictorNames( ...
        ~ismember( ...
        requestedPredictorNames, ...
        availablePredictorNames));

    if ~isempty(invalidPredictorNames)
        error( ...
            "plotFullReducedModelPETHs:InvalidPredictor", ...
            ["The following requested predictors are unavailable, " ...
            "lack a fitted reduced model, or are not discrete-lagged: %s"], ...
            strjoin(invalidPredictorNames,", "));
    end

    predictorNames = requestedPredictorNames;
end

if isempty(predictorNames)
    error( ...
        "plotFullReducedModelPETHs:NoPredictors", ...
        ["No discrete-lagged predictors with fitted reduced models " ...
        "were found."]);
end

nPredictors = numel(predictorNames);

resultPredictorIndex = nan(1,nPredictors);

for predictorIndex = 1:nPredictors

    idx = find( ...
        resultPredictorNames == predictorNames(predictorIndex), ...
        1);

    if isempty(idx)
        error( ...
            "plotFullReducedModelPETHs:UnknownPredictor", ...
            "Predictor '%s' was not found in results.predictor.", ...
            predictorNames(predictorIndex));
    end

    if isempty(results.predictor(idx).otherModel)
        error( ...
            "plotFullReducedModelPETHs:MissingReducedModel", ...
            "The reduced model for predictor '%s' is empty.", ...
            predictorNames(predictorIndex));
    end

    resultPredictorIndex(predictorIndex) = idx;
end

%% Compute empirical and full-model PETHs

empiricalPETH = computeConditionPETHs( ...
    P, ...
    F, ...
    predictorInfo);

FpredFull = predictPoissonGLM( ...
    results.full.model, ...
    P);

fullPredictedPETH = computeConditionPETHs( ...
    P, ...
    FpredFull, ...
    predictorInfo);

%% Compute reduced-model predictions and PETHs

reducedPredictedPETH = cell(1,nPredictors);

for predictorIndex = 1:nPredictors

    resultIndex = resultPredictorIndex(predictorIndex);
    predictorResult = results.predictor(resultIndex);

    if ~isfield(predictorResult,"otherColumns") || ...
            isempty(predictorResult.otherColumns)

        error( ...
            "plotFullReducedModelPETHs:MissingOtherColumns", ...
            ["results.predictor(%d).otherColumns is required to " ...
            "construct predictions from the reduced model."], ...
            resultIndex);
    end

    otherColumns = predictorResult.otherColumns;

    if any(otherColumns < 1) || any(otherColumns > size(P,2))
        error( ...
            "plotFullReducedModelPETHs:InvalidOtherColumns", ...
            "Invalid design-matrix column indices for predictor '%s'.", ...
            predictorNames(predictorIndex));
    end

    Pred = P(:,otherColumns);

    FpredReduced = predictPoissonGLM( ...
        predictorResult.otherModel, ...
        Pred);

    % Use the full P and predictorInfo to define trial conditions and
    % alignment. Only the predicted activity differs.
    reducedPredictedPETH{predictorIndex} = ...
        computeConditionPETHs( ...
        P, ...
        FpredReduced, ...
        predictorInfo);
end

%% Match PETH groups

empiricalNames = string({empiricalPETH.name});
fullNames = string({fullPredictedPETH.name});

empiricalIndex = nan(1,nPredictors);
fullIndex = nan(1,nPredictors);
reducedIndex = nan(1,nPredictors);

for predictorIndex = 1:nPredictors

    predictorName = predictorNames(predictorIndex);

    idx = find(empiricalNames == predictorName,1);
    if ~isempty(idx)
        empiricalIndex(predictorIndex) = idx;
    end

    idx = find(fullNames == predictorName,1);
    if ~isempty(idx)
        fullIndex(predictorIndex) = idx;
    end

    reducedNames = string( ...
        {reducedPredictedPETH{predictorIndex}.name});

    idx = find(reducedNames == predictorName,1);
    if ~isempty(idx)
        reducedIndex(predictorIndex) = idx;
    end
end

%% Make one figure per neuron

for neuronIndex = neuronIndices(:).'

    figureWidth = max(520,130*nPredictors);

    fig = figure( ...
        "WindowStyle","normal", ...
        "Name",sprintf( ...
        "Neuron %d: empirical, full, and reduced PETHs", ...
        neuronIndex), ...
        "Color","w", ...
        "Position",[100 200 figureWidth 350]);

    layout = tiledlayout( ...
        fig, ...
        3, ...
        nPredictors, ...
        "TileSpacing","compact", ...
        "Padding","compact");

    axesHandles = gobjects(3,nPredictors);

    for predictorIndex = 1:nPredictors

        predictorName = predictorNames(predictorIndex);

        %% First row: empirical activity

        axEmpirical = nexttile( ...
            layout, ...
            predictorIndex);

        axesHandles(1,predictorIndex) = axEmpirical;

        plotOneModelPanel( ...
            axEmpirical, ...
            empiricalPETH, ...
            empiricalPETH, ...
            empiricalIndex(predictorIndex), ...
            empiricalIndex(predictorIndex), ...
            predictorName, ...
            neuronIndex, ...
            "empirical", ...
            options);

        title( ...
            axEmpirical, ...
            replace(predictorName,"_"," "));

        %% Second row: full-model prediction

        axFull = nexttile( ...
            layout, ...
            nPredictors+predictorIndex);

        axesHandles(2,predictorIndex) = axFull;

        plotOneModelPanel( ...
            axFull, ...
            empiricalPETH, ...
            fullPredictedPETH, ...
            empiricalIndex(predictorIndex), ...
            fullIndex(predictorIndex), ...
            predictorName, ...
            neuronIndex, ...
            "full", ...
            options);

        %% Third row: reduced-model prediction

        axReduced = nexttile( ...
            layout, ...
            2*nPredictors+predictorIndex);

        axesHandles(3,predictorIndex) = axReduced;

        thisReducedPETH = ...
            reducedPredictedPETH{predictorIndex};

        plotOneModelPanel( ...
            axReduced, ...
            empiricalPETH, ...
            thisReducedPETH, ...
            empiricalIndex(predictorIndex), ...
            reducedIndex(predictorIndex), ...
            predictorName, ...
            neuronIndex, ...
            "reduced", ...
            options);

        %Add unique cross-validated deviance explained
        resultIndex = resultPredictorIndex(predictorIndex);
        cvDE = results.deUnique(neuronIndex,resultIndex);
        if isfinite(cvDE)
            if 100*cvDE<=1
                txtCol = [.6 .6 .6];
            else
                txtCol = [0 0 0];
            end
            text( ...
                axReduced, ...
                0.98, ...
                0.95, ...
                sprintf("cvDE=%.1f%%",100*cvDE), ...
                "Units","normalized", ...
                "HorizontalAlignment","right", ...
                "VerticalAlignment","top", ...
                "Color",txtCol, ...
                "Margin",2, ...
                "FontSize",8, ...
                "Clipping","on");
        else
            text( ...
                axReduced, ...
                0.98, ...
                0.95, ...
                "cvDE=n/a", ...
                "Units","normalized", ...
                "HorizontalAlignment","right", ...
                "VerticalAlignment","top", ...
                "Color",[0.6 0.6 0.6], ...
                "Margin",2, ...
                "FontSize",8, ...
                "Clipping","on");
        end

    end

    %% Axis labels

    ylabel( ...
        axesHandles(1,1), ...
        sprintf("Empirical\nactivity"));

    ylabel( ...
        axesHandles(2,1), ...
        sprintf("Full model\nactivity"));

    ylabel( ...
        axesHandles(3,1), ...
        sprintf("Reduced model\nactivity"));

    % Hide redundant y-axis tick labels.
    for rowIndex = 1:3
        for predictorIndex = 2:nPredictors
            axesHandles(rowIndex,predictorIndex).YTickLabel = [];
        end
    end

    % Show x-axis labels only on the bottom row.
    for predictorIndex = 1:nPredictors

        axesHandles(1,predictorIndex).XTickLabel = [];
        axesHandles(2,predictorIndex).XTickLabel = [];

        xlabel( ...
            axesHandles(3,predictorIndex), ...
            "Lag (s)");
    end

    %% Match y-limits

    if options.shareYAcrossRows

        validAxes = axesHandles(isgraphics(axesHandles,"axes"));

        linkaxes(validAxes,"y");

        allLimits = vertcat(validAxes.YLim);

        commonYLim = [ ...
            min(allLimits(:,1)), ...
            max(allLimits(:,2))];

        if commonYLim(1) >= 0
            commonYLim(1) = 0;
        end

        set(validAxes,"YLim",commonYLim);

    else

        % Share y-limits across predictors within each row.
        for rowIndex = 1:3
            linkaxes(axesHandles(rowIndex,:),"y");
        end
    end

    title( ...
        layout, ...
        sprintf("Neuron %d: data, full model, and predictor-omitted models", ...
        neuronIndex));
end

end

function lineHandles = plotOneModelPanel( ...
    ax, empiricalPETH, predictedPETH, ...
    empiricalGroupIndex, predictedGroupIndex, ...
    predictorName, neuronIndex, panelType, options)
%PLOTONEMODELPANEL Plot an empirical or model-predicted PETH panel.
%
%   panelType = "empirical"
%       Empirical mean as a solid line with empirical SEM shading.
%
%   panelType = "full"
%       Full-model prediction as a dashed line.
%
%   panelType = "reduced"
%       Predictor-omitted model prediction as a dashed line.

arguments
    ax (1,1) matlab.graphics.axis.Axes
    empiricalPETH struct
    predictedPETH struct
    empiricalGroupIndex (1,1) double
    predictedGroupIndex (1,1) double
    predictorName (1,1) string
    neuronIndex (1,1) double {mustBeInteger,mustBePositive}
    panelType (1,1) string ...
        {mustBeMember(panelType,["empirical","full","reduced"])}
    options struct
end

hold(ax,"on");

%% Handle unavailable PETH groups

if isnan(empiricalGroupIndex) || isnan(predictedGroupIndex)

    text( ...
        ax, ...
        0.5, ...
        0.5, ...
        "Matching PETH not found", ...
        "Units","normalized", ...
        "HorizontalAlignment","center", ...
        "VerticalAlignment","middle");

    box(ax,"off");

    lineHandles = gobjects(0);
    return
end

empiricalGroup = empiricalPETH(empiricalGroupIndex);
predictedGroup = predictedPETH(predictedGroupIndex);

%% Verify neuron index

nEmpiricalNeurons = size(empiricalGroup.mean,1);
nPredictedNeurons = size(predictedGroup.mean,1);

if neuronIndex > nEmpiricalNeurons
    error( ...
        "plotFullReducedModelPETHs:EmpiricalNeuronIndex", ...
        ["Neuron %d exceeds the %d neurons available in the " ...
        "empirical PETH for '%s'."], ...
        neuronIndex, ...
        nEmpiricalNeurons, ...
        predictorName);
end

if neuronIndex > nPredictedNeurons
    error( ...
        "plotFullReducedModelPETHs:PredictedNeuronIndex", ...
        ["Neuron %d exceeds the %d neurons available in the " ...
        "predicted PETH for '%s'."], ...
        neuronIndex, ...
        nPredictedNeurons, ...
        predictorName);
end

%% Verify empirical and predicted PETH alignment

empiricalLags = empiricalGroup.lags(:).';
predictedLags = predictedGroup.lags(:).';

empiricalConditions = empiricalGroup.conditionValues(:).';
predictedConditions = predictedGroup.conditionValues(:).';

if numel(empiricalLags) ~= numel(predictedLags) || ...
        any(abs(empiricalLags-predictedLags) > 1e-10)

    error( ...
        "plotFullReducedModelPETHs:LagMismatch", ...
        "Empirical and predicted lags do not match for '%s'.", ...
        predictorName);
end

if numel(empiricalConditions) ~= numel(predictedConditions) || ...
        any(abs(empiricalConditions-predictedConditions) > 1e-10)

    error( ...
        "plotFullReducedModelPETHs:ConditionMismatch", ...
        ["Empirical and predicted condition values do not match " ...
        "for '%s'."], ...
        predictorName);
end

lags = empiricalLags;
conditionValues = empiricalConditions;

nLags = numel(lags);
nConditions = numel(conditionValues);

%% Extract this neuron's PETHs

empiricalMean = reshape( ...
    empiricalGroup.mean(neuronIndex,:,:), ...
    nLags, ...
    nConditions);

empiricalSEM = reshape( ...
    empiricalGroup.sem(neuronIndex,:,:), ...
    nLags, ...
    nConditions);

predictedMean = reshape( ...
    predictedGroup.mean(neuronIndex,:,:), ...
    nLags, ...
    nConditions);

%% Optional baseline subtraction

if options.baselineSubtract

    baselineMask = ...
        lags >= options.baselineWindow(1) & ...
        lags <= options.baselineWindow(2);

    if ~any(baselineMask)
        error( ...
            "plotFullReducedModelPETHs:InvalidBaselineWindow", ...
            "No bins for '%s' fall within [%.3g, %.3g] seconds.", ...
            predictorName, ...
            options.baselineWindow(1), ...
            options.baselineWindow(2));
    end

    empiricalBaseline = mean( ...
        empiricalMean(baselineMask,:), ...
        1, ...
        "omitnan");

    predictedBaseline = mean( ...
        predictedMean(baselineMask,:), ...
        1, ...
        "omitnan");

    empiricalMean = empiricalMean-empiricalBaseline;
    predictedMean = predictedMean-predictedBaseline;
end

%% Select line source and style

switch panelType

    case "empirical"
        lineValues = empiricalMean;
        lineStyle = "-";
        lineWidth = options.empiricalLineWidth;

    case "full"
        lineValues = predictedMean;
        lineStyle = "-";
        lineWidth = options.fullLineWidth;

    case "reduced"
        lineValues = predictedMean;
        lineStyle = "--";
        lineWidth = options.reducedLineWidth;
end

%% Set condition colors

colorOrder = ax.ColorOrder;

contrastSaturation = [0.05,0.25,0.65,1.00];
contrastValue = [0.7,0.9,0.9,0.7];

leftHue = 0.00;
rightHue = 0.58;

leftContrastColors = hsv2rgb([ ...
    repmat(leftHue,4,1), ...
    contrastSaturation(:), ...
    contrastValue(:)]);

rightContrastColors = hsv2rgb([ ...
    repmat(rightHue,4,1), ...
    contrastSaturation(:), ...
    contrastValue(:)]);

lineHandles = gobjects(nConditions,1);

%% Plot each condition

for conditionIndex = 1:nConditions

    switch predictorName

        case "stimulusLeft"
            thisColor = leftContrastColors(conditionIndex,:);

        case "stimulusRight"
            thisColor = rightContrastColors(conditionIndex,:);

        otherwise
            colorIndex = ...
                mod(conditionIndex-1,size(colorOrder,1))+1;

            thisColor = colorOrder(colorIndex,:);
    end

    thisEmpiricalMean = ...
        empiricalMean(:,conditionIndex);

    thisEmpiricalSEM = ...
        empiricalSEM(:,conditionIndex);

    thisLineValues = ...
        lineValues(:,conditionIndex);

    %% Empirical SEM shading: first row only

    if panelType == "empirical" && options.showSEM

        lowerBound = ...
            thisEmpiricalMean-thisEmpiricalSEM;

        upperBound = ...
            thisEmpiricalMean+thisEmpiricalSEM;

        fill( ...
            ax, ...
            [lags,fliplr(lags)], ...
            [lowerBound.',fliplr(upperBound.')], ...
            thisColor, ...
            "EdgeColor","none", ...
            "FaceAlpha",options.semAlpha, ...
            "HandleVisibility","off");
    end

    %% Empirical or predicted mean

    lineHandles(conditionIndex) = plot( ...
        ax, ...
        lags, ...
        thisLineValues, ...
        lineStyle, ...
        "Color",thisColor, ...
        "LineWidth",lineWidth);
end

%% Format axes

xline( ...
    ax, ...
    0, ...
    "k--", ...
    "LineWidth",1, ...
    "HandleVisibility","off");

xlim(ax,[min(lags),max(lags)]);

box(ax,"off");

%% Optional condition legend

% Show the condition legend only in the empirical row to avoid repeating
% the same legend three times.
if options.showLegend && panelType == "empirical"

    conditionLabels = makeConditionLabels( ...
        predictorName, ...
        conditionValues);

    conditionLegend = legend( ...
        ax, ...
        flipud(lineHandles), ...
        flip(conditionLabels), ...
        "Location","best", ...
        "Box","off");

    conditionLegend.AutoUpdate = "off";
end

end