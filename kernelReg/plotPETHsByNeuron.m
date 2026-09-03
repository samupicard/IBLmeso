function plotPETHsByNeuron(peth, neuronIndices, options)
%PLOTPETHSBYNEURON Plot empirical PETHs for selected neurons.
%
% One figure is created per neuron, with panels for:
%   1. stimulusLeft
%   2. stimulusRight
%   3. choice
%   4. feedback
%
% Expected PETH structure:
%   peth(i).name
%   peth(i).lags                 1 × nLags
%   peth(i).mean                 neurons × lags × conditions
%   peth(i).sem                  neurons × lags × conditions
%   peth(i).conditionValues      1 × conditions
%
% Example:
%   plotPETHsByNeuron( ...
%       empiricalPETH, ...
%       [12 48 103], ...
%       showSEM = true);

arguments
    peth struct
    neuronIndices (:,1) double {mustBeInteger, mustBePositive}

    options.showSEM (1,1) logical = true
    options.baselineSubtract (1,1) logical = false
    options.baselineWindow (1,2) double = [-1 -0.5]
    options.lineWidth (1,1) double = 1.5
end

predictorNames = [ ...
    "stimulusLeft", ...
    "stimulusRight", ...
    "choice", ...
    "feedback"];

panelTitles = [ ...
    "Stimulus left", ...
    "Stimulus right", ...
    "Choice", ...
    "Feedback"];

%% Find requested predictor groups

pethNames = string({peth.name});

pethIndex = nan(size(predictorNames));

for ii = 1:numel(predictorNames)
    idx = find(pethNames == predictorNames(ii), 1);

    if ~isempty(idx)
        pethIndex(ii) = idx;
    end
end

%% Check neuron indices

availableIndex = find(~isnan(pethIndex), 1);

if isempty(availableIndex)
    error("None of the requested PETH groups were found.");
end

examplePETH = peth(pethIndex(availableIndex));
nNeurons = size(examplePETH.mean, 1);

if any(neuronIndices > nNeurons)
    error( ...
        "Requested neuron index exceeds the available %d neurons.", ...
        nNeurons);
end

%% Create one figure per neuron

for neuronIndex = neuronIndices(:).'
    fig = figure( ...
        "Name", sprintf("Neuron %d empirical PETHs", neuronIndex), ...
        "Color", "w");

    layout = tiledlayout( ...
        fig, ...
        1, 4, ...
        "TileSpacing", "compact", ...
        "Padding", "compact");

    axesHandles = gobjects(1, numel(predictorNames));

    for panelIndex = 1:numel(predictorNames)
        ax = nexttile(layout);
        axesHandles(panelIndex) = ax;

        hold(ax, "on");

        thisPETHIndex = pethIndex(panelIndex);

        if isnan(thisPETHIndex)
            title(ax, panelTitles(panelIndex));

            text( ...
                ax, ...
                0.5, ...
                0.5, ...
                "Predictor not found", ...
                "Units", "normalized", ...
                "HorizontalAlignment", "center");

            continue
        end

        thisPETH = peth(thisPETHIndex);

        lags = thisPETH.lags(:).';
        conditionValues = thisPETH.conditionValues(:).';

        meanResponse = reshape( ...
            thisPETH.mean(neuronIndex, :, :), ...
            numel(lags), ...
            numel(conditionValues));

        semResponse = reshape( ...
            thisPETH.sem(neuronIndex, :, :), ...
            numel(lags), ...
            numel(conditionValues));

        %% Optional baseline subtraction

        if options.baselineSubtract
            baselineMask = ...
                lags >= options.baselineWindow(1) & ...
                lags <= options.baselineWindow(2);

            if ~any(baselineMask)
                error( ...
                    "No lag bins for %s fall within [%.3g, %.3g].", ...
                    predictorNames(panelIndex), ...
                    options.baselineWindow(1), ...
                    options.baselineWindow(2));
            end

            baseline = mean( ...
                meanResponse(baselineMask, :), ...
                1, ...
                "omitnan");

            meanResponse = meanResponse - baseline;
        end

        %% Plot conditions

        lineHandles = gobjects(numel(conditionValues), 1);
        colorOrder = ax.ColorOrder;

        for conditionIndex = 1:numel(conditionValues)
            thisMean = meanResponse(:, conditionIndex);
            thisSEM = semResponse(:, conditionIndex);

            thisColor = colorOrder( ...
                mod(conditionIndex - 1, size(colorOrder, 1)) + 1, :);

            if options.showSEM
                fill( ...
                    ax, ...
                    [lags, fliplr(lags)], ...
                    [ ...
                    (thisMean - thisSEM).', ...
                    fliplr((thisMean + thisSEM).') ...
                    ], ...
                    thisColor, ...
                    "EdgeColor", "none", ...
                    "FaceAlpha", 0.15, ...
                    "HandleVisibility", "off");
            end

            lineHandles(conditionIndex) = plot( ...
                ax, ...
                lags, ...
                thisMean, ...
                "Color", thisColor, ...
                "LineWidth", options.lineWidth);
        end

        xline( ...
            ax, ...
            0, ...
            "k--", ...
            "LineWidth", 1, ...
            "HandleVisibility", "off");

        xlabel(ax, "Lag (s)");
        title(ax, panelTitles(panelIndex));

        conditionLabels = makeConditionLabels( ...
            predictorNames(panelIndex), ...
            conditionValues);

        legend( ...
            ax, ...
            lineHandles, ...
            conditionLabels, ...
            "Location", "best", ...
            "Box", "off");

        box(ax, "off");
    end

    %% Share the y-axis within this neuron

    validAxes = axesHandles(isgraphics(axesHandles, "axes"));
    linkaxes(validAxes, "y");

    ylabel(validAxes(1), "Neural activity");

    for axIndex = 2:numel(validAxes)
        ylabel(validAxes(axIndex), "");
        validAxes(axIndex).YTickLabel = [];
    end

    title(layout, sprintf("Neuron %d: empirical PETHs", neuronIndex));
end

end


function labels = makeConditionLabels(predictorName, conditionValues)
%MAKECONDITIONLABELS Generate readable condition labels.

switch predictorName
    case {"stimulusLeft", "stimulusRight"}
        labels = compose("Contrast %g", conditionValues);

    case "choice"
        labels = strings(size(conditionValues));

        for ii = 1:numel(conditionValues)
            switch conditionValues(ii)
                case -1
                    labels(ii) = "Left choice";
                case 1
                    labels(ii) = "Right choice";
                otherwise
                    labels(ii) = sprintf("Choice %g", conditionValues(ii));
            end
        end

    case "feedback"
        labels = strings(size(conditionValues));

        for ii = 1:numel(conditionValues)
            switch conditionValues(ii)
                case -1
                    labels(ii) = "Incorrect";
                case 1
                    labels(ii) = "Correct";
                otherwise
                    labels(ii) = sprintf("Feedback %g", conditionValues(ii));
            end
        end

    otherwise
        labels = compose("%g", conditionValues);
end

end