function plotEmpiricalPredictedPETHs( ...
    empiricalPETH, predictedPETH, neuronIndices, options)
%PLOTEMPIRICALPREDICTEDPETHS Compare observed and predicted PETHs.
%
% One figure is produced per neuron, containing four panels:
%   stimulusLeft | stimulusRight | choice | feedback
%
% For each condition:
%   Empirical mean   = solid line
%   Empirical SEM    = shaded region
%   Predicted mean   = dashed line
%
% Inputs
%   empiricalPETH
%       Output from computeConditionPETHs(P, F, predictorInfo).
%
%   predictedPETH
%       Output from computeConditionPETHs(P, Fpred, predictorInfo).
%
%   neuronIndices
%       Indices of neurons to plot.
%
% Options
%   showSEM
%       Display empirical SEM shading.
%
%   baselineSubtract
%       Baseline-subtract both empirical and predicted PETHs.
%
%   baselineWindow
%       Baseline interval in seconds.
%
%   empiricalLineWidth
%       Width of empirical solid lines.
%
%   predictedLineWidth
%       Width of predicted dashed lines.

arguments
    empiricalPETH struct
    predictedPETH struct
    neuronIndices (:,1) double {mustBeInteger, mustBePositive}

    options.showSEM (1,1) logical = true
    options.showLegend (1,1) logical = false
    options.baselineSubtract (1,1) logical = false
    options.baselineWindow (1,2) double = [-1 -0.5]

    options.empiricalLineWidth (1,1) double = 1.75
    options.predictedLineWidth (1,1) double = 1.75
    options.semAlpha (1,1) double = 0.15
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

%% Match predictor groups

empiricalNames = string({empiricalPETH.name});
predictedNames = string({predictedPETH.name});

empiricalIndex = nan(size(predictorNames));
predictedIndex = nan(size(predictorNames));

for ii = 1:numel(predictorNames)
    idx = find(empiricalNames == predictorNames(ii), 1);

    if ~isempty(idx)
        empiricalIndex(ii) = idx;
    end

    idx = find(predictedNames == predictorNames(ii), 1);

    if ~isempty(idx)
        predictedIndex(ii) = idx;
    end
end

availablePanel = find( ...
    ~isnan(empiricalIndex) & ~isnan(predictedIndex), ...
    1);

if isempty(availablePanel)
    error( ...
        "No matching predictor groups were found in the empirical " + ...
        "and predicted PETH structures.");
end

nNeuronsEmpirical = size( ...
    empiricalPETH(empiricalIndex(availablePanel)).mean, 1);

nNeuronsPredicted = size( ...
    predictedPETH(predictedIndex(availablePanel)).mean, 1);

if nNeuronsEmpirical ~= nNeuronsPredicted
    error( ...
        "Empirical and predicted PETHs have different neuron counts.");
end

if any(neuronIndices > nNeuronsEmpirical)
    error( ...
        "A requested neuron index exceeds the available %d neurons.", ...
        nNeuronsEmpirical);
end

%% Make one figure per neuron

for neuronIndex = neuronIndices(:).'
    fig = figure( ...
        "WindowStyle", "normal", ...
        "Name", sprintf( ...
        "Neuron %d: empirical and predicted PETHs", ...
        neuronIndex), ...
        "Color", "w",...
        "Position",[160,100,500,300]);

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

        empiricalGroupIndex = empiricalIndex(panelIndex);
        predictedGroupIndex = predictedIndex(panelIndex);

        if isnan(empiricalGroupIndex) || isnan(predictedGroupIndex)
            title(ax, panelTitles(panelIndex));

            text( ...
                ax, ...
                0.5, ...
                0.5, ...
                "Matching PETH not found", ...
                "Units", "normalized", ...
                "HorizontalAlignment", "center");

            continue
        end

        empiricalGroup = empiricalPETH(empiricalGroupIndex);
        predictedGroup = predictedPETH(predictedGroupIndex);

        %% Verify that empirical and predicted PETHs align

        empiricalLags = empiricalGroup.lags(:).';
        predictedLags = predictedGroup.lags(:).';

        empiricalConditions = empiricalGroup.conditionValues(:).';
        predictedConditions = predictedGroup.conditionValues(:).';

        if numel(empiricalLags) ~= numel(predictedLags) || ...
                any(abs(empiricalLags - predictedLags) > 1e-10)

            error( ...
                "Empirical and predicted lags do not match for '%s'.", ...
                predictorNames(panelIndex));
        end

        if numel(empiricalConditions) ~= numel(predictedConditions) || ...
                any(abs( ...
                empiricalConditions - predictedConditions) > 1e-10)

            error( ...
                ["Empirical and predicted condition values do not " ...
                "match for '%s'."], ...
                predictorNames(panelIndex));
        end

        lags = empiricalLags;
        conditionValues = empiricalConditions;

        nLags = numel(lags);
        nConditions = numel(conditionValues);

        empiricalMean = reshape( ...
            empiricalGroup.mean(neuronIndex, :, :), ...
            nLags, ...
            nConditions);

        empiricalSEM = reshape( ...
            empiricalGroup.sem(neuronIndex, :, :), ...
            nLags, ...
            nConditions);

        predictedMean = reshape( ...
            predictedGroup.mean(neuronIndex, :, :), ...
            nLags, ...
            nConditions);

        %% Optional baseline subtraction

        if options.baselineSubtract
            baselineMask = ...
                lags >= options.baselineWindow(1) & ...
                lags <= options.baselineWindow(2);

            if ~any(baselineMask)
                error( ...
                    "No bins for %s fall within [%.3g, %.3g].", ...
                    predictorNames(panelIndex), ...
                    options.baselineWindow(1), ...
                    options.baselineWindow(2));
            end

            % Subtract each curve's own baseline.
            empiricalBaseline = mean( ...
                empiricalMean(baselineMask, :), ...
                1, ...
                "omitnan");

            predictedBaseline = mean( ...
                predictedMean(baselineMask, :), ...
                1, ...
                "omitnan");

            empiricalMean = empiricalMean - empiricalBaseline;
            predictedMean = predictedMean - predictedBaseline;
        end

        %% Plot conditions

        empiricalHandles = gobjects(nConditions,1);

        colorOrder = ax.ColorOrder;

        % Default MATLAB colors
        blue = colorOrder(1,:);   % first color
        red  = colorOrder(2,:);   % second color

        switch predictorNames(panelIndex)

            case "stimulusLeft"
                baseColor = red;

            case "stimulusRight"
                baseColor = blue;

            otherwise
                baseColor = [];
        end

        for conditionIndex = 1:nConditions
            if isempty(baseColor)

                % Original color scheme for choice/feedback
                thisColor = colorOrder( ...
                    mod(conditionIndex-1,size(colorOrder,1))+1,:);

            else

                % Strongest contrast = darkest.
                % Blend from baseColor toward white.
                % conditionIndex = n  -> baseColor
                % conditionIndex = 1  -> very light version

                blend = (nConditions-conditionIndex)/max(nConditions-1,1);

                % don't go all the way to white
                blend = 0.8*blend;

                thisColor = (1-blend)*baseColor + blend*[1 1 1];

            end
            thisEmpiricalMean = empiricalMean(:, conditionIndex);
            thisEmpiricalSEM = empiricalSEM(:, conditionIndex);
            thisPredictedMean = predictedMean(:, conditionIndex);

            % Empirical SEM region.
            if options.showSEM
                fill( ...
                    ax, ...
                    [lags, fliplr(lags)], ...
                    [ ...
                    (thisEmpiricalMean - thisEmpiricalSEM).', ...
                    fliplr( ...
                    (thisEmpiricalMean + ...
                    thisEmpiricalSEM).') ...
                    ], ...
                    thisColor, ...
                    "EdgeColor", "none", ...
                    "FaceAlpha", options.semAlpha, ...
                    "HandleVisibility", "off");
            end

            % Empirical mean.
            empiricalHandles(conditionIndex) = plot( ...
                ax, ...
                lags, ...
                thisEmpiricalMean, ...
                "-", ...
                "Color", thisColor, ...
                "LineWidth", options.empiricalLineWidth);

            % Predicted mean.
            plot( ...
                ax, ...
                lags, ...
                thisPredictedMean, ...
                "--", ...
                "Color", thisColor, ...
                "LineWidth", options.predictedLineWidth, ...
                "HandleVisibility", "off");
        end

        xline( ...
            ax, ...
            0, ...
            "k--", ...
            "LineWidth", 1, ...
            "HandleVisibility", "off");

        xlim(ax, [min(lags), max(lags)]);

        xlabel(ax, "Lag (s)");
        title(ax, panelTitles(panelIndex));
        box(ax, "off");

        if options.showLegend
        conditionLabels = makeConditionLabels( ...
            predictorNames(panelIndex), ...
            conditionValues);

        % Color legend: one entry per condition.
        conditionLegend = legend( ...
            ax, ...
            flipud(empiricalHandles), ...
            flip(conditionLabels), ...
            "Location", "best", ...
            "Box", "off");

        conditionLegend.AutoUpdate = "off";
        end
    end

    %% Share y-axis across all panels for this neuron

    validAxes = axesHandles(isgraphics(axesHandles, "axes"));
    linkaxes(validAxes, "y");

    ylims = get(gca,'ylim');
    ylim([0,ylims(end)]);

    ylabel(validAxes(1), "Neural activity");

    for axIndex = 2:numel(validAxes)
        ylabel(validAxes(axIndex), "");
        validAxes(axIndex).YTickLabel = [];
    end

    title(layout, sprintf( ...
        "Neuron %d: empirical (solid) and predicted (dashed)", ...
        neuronIndex));
end

end


function labels = makeConditionLabels(predictorName, conditionValues)
%MAKECONDITIONLABELS Produce readable condition labels.

switch predictorName
    case {"stimulusLeft", "stimulusRight"}
        labels = compose("%g", conditionValues);

    case "choice"
        labels = strings(size(conditionValues));

        for ii = 1:numel(conditionValues)
            switch conditionValues(ii)
                case -1
                    labels(ii) = "Right";
                case 1
                    labels(ii) = "Left";
                otherwise
                    labels(ii) = sprintf( ...
                        "%g", conditionValues(ii));
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
                    labels(ii) = sprintf( ...
                        "Feedback %g", conditionValues(ii));
            end
        end

    otherwise
        labels = compose("%g", conditionValues);
end

end