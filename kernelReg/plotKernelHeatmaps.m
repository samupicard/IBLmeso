function plotKernelHeatmaps(model, predictorInfo, options)
%PLOTKERNELHEATMAPS Plot GLM kernels and non-kernel predictor weights.
%
% model.K:
%   predictors × neurons, excluding the intercept.
%
% predictorInfo:
%   Table with variables:
%       column
%       name
%       type
%       lagSeconds
%
% Lagged predictors are plotted as neuron × time heatmaps.
%
% Predictors with type "epoch" or "continuous" are plotted as narrow,
% single-column heatmaps using the same neuron ordering and colour scale.
%
% SORT MODES
% ----------
% "cluster"
%   Orders neurons by hierarchical clustering of their coefficient profiles
%   across all displayed predictors. Kernel groups are scaled so that
%   predictors with many lag bins do not automatically dominate the
%   clustering.
%
% "peakLag"
%   Orders neurons by the lag of the largest absolute coefficient in the
%   first lagged predictor.
%
% "kernelStrength"
%   Orders neurons by the L2 norm of the first lagged predictor.
%
% "none"
%   Preserves the original neuron order.
%
% Example:
%   plotKernelHeatmaps(model, predictorInfo, ...
%       minimumDevianceExplained=0.05, ...
%       sortMode="cluster");

arguments
    model struct
    predictorInfo table

    options.sortMode (1,1) string ...
        {mustBeMember(options.sortMode, ...
        ["cluster","peakLag","kernelStrength","none"])} = "peakLag"

    options.minimumDevianceExplained (1,1) double = -Inf

    options.standardize (1,1) logical = false
    options.P double = double.empty

    options.showNeuronNumbers (1,1) logical = false

    options.colormapName double = brewermap(256,'*RdBu')
    options.percentileLimit (1,1) double ...
        {mustBeGreaterThan(options.percentileLimit,0), ...
         mustBeLessThanOrEqual(options.percentileLimit,100)} = 99

    % Relative tile widths. Lagged predictors span several tile columns,
    % whereas scalar epoch/continuous predictors span one.
    options.kernelTileWidth (1,1) double ...
        {mustBeInteger,mustBePositive} = 4

    options.scalarTileWidth (1,1) double ...
        {mustBeInteger,mustBePositive} = 1

    options.figureHeight (1,1) double ...
        {mustBePositive} = 650

    options.kernelPanelWidth (1,1) double ...
        {mustBePositive} = 60

    options.scalarPanelWidth (1,1) double ...
        {mustBePositive} = 30

    options.minimumFigureWidth (1,1) double ...
        {mustBePositive} = 400
end

requiredVariables = ["column","name","type","lagSeconds"];

if ~all(ismember(requiredVariables,string(predictorInfo.Properties.VariableNames)))
    error( ...
        "plotKernelHeatmaps:InvalidPredictorInfo", ...
        "predictorInfo must contain variables: %s.", ...
        strjoin(requiredVariables,", "));
end

assert( ...
    size(model.K,1) == height(predictorInfo), ...
    "model.K rows must correspond directly to predictorInfo rows.");

%% Select well-fit neurons

if ~isfield(model,"devianceExplained")
    error( ...
        "plotKernelHeatmaps:MissingDevianceExplained", ...
        "model must contain the field model.devianceExplained.");
end

devianceExplained = model.devianceExplained(:);

if numel(devianceExplained) ~= size(model.K,2)
    error( ...
        "plotKernelHeatmaps:NeuronCountMismatch", ...
        "model.devianceExplained must contain one value per neuron.");
end

neuronMask = ...
    isfinite(devianceExplained) & ...
    devianceExplained >= options.minimumDevianceExplained;

if ~any(neuronMask)
    error( ...
        "plotKernelHeatmaps:NoNeurons", ...
        "No neurons satisfy devianceExplained >= %.3g.", ...
        options.minimumDevianceExplained);
end

selectedNeuronIndices = find(neuronMask);
K = double(model.K(:,neuronMask));

nSelectedNeurons = size(K,2);

%% Optionally standardize coefficients

if options.standardize

    if isempty(options.P)
        error( ...
            "plotKernelHeatmaps:MissingPredictorMatrix", ...
            "Supply options.P when standardize=true.");
    end

    if size(options.P,2) ~= size(K,1)
        error( ...
            "plotKernelHeatmaps:PredictorCountMismatch", ...
            "options.P must contain one column per row of model.K.");
    end

    predictorScale = std(options.P,0,1,"omitnan");
    predictorScale(~isfinite(predictorScale)) = 0;

    K = K .* predictorScale(:);
end

%% Identify lagged and scalar predictors

predictorTypes = string(predictorInfo.type);
predictorNames = string(predictorInfo.name);

isLagged = ...
    predictorTypes == "discrete-lagged" & ...
    isfinite(predictorInfo.lagSeconds);

isScalar = ...
    predictorTypes == "epoch" | ...
    predictorTypes == "continuous";

laggedInfo = predictorInfo(isLagged,:);
kernelNames = unique(string(laggedInfo.name),"stable");

scalarRows = find(isScalar);

nKernelGroups = numel(kernelNames);
nScalarPredictors = numel(scalarRows);
nPanels = nKernelGroups+nScalarPredictors;

if nPanels == 0
    error( ...
        "plotKernelHeatmaps:NoPredictors", ...
        "No discrete-lagged, epoch, or continuous predictors were found.");
end

%% Construct kernel groups

kernelRows = cell(nKernelGroups,1);
kernelLags = cell(nKernelGroups,1);
kernelValues = cell(nKernelGroups,1);

for kernelIndex = 1:nKernelGroups

    groupName = kernelNames(kernelIndex);
    groupMask = string(laggedInfo.name) == groupName;

    rows = laggedInfo.column(groupMask);
    lags = laggedInfo.lagSeconds(groupMask);

    [lags,lagOrder] = sort(lags);
    rows = rows(lagOrder);

    kernelRows{kernelIndex} = rows(:);
    kernelLags{kernelIndex} = lags(:).';

    % Neurons × lag bins
    kernelValues{kernelIndex} = K(rows,:).';
end

%% Construct scalar predictor values and labels

scalarValues = cell(nScalarPredictors,1);
scalarLabels = strings(nScalarPredictors,1);

for scalarIndex = 1:nScalarPredictors

    infoRow = scalarRows(scalarIndex);
    coefficientRow = predictorInfo.column(infoRow);

    scalarValues{scalarIndex} = K(coefficientRow,:).';

    scalarLabels(scalarIndex) = string(predictorNames{infoRow}(1:2));

    % Add the type only when duplicate names would otherwise be ambiguous.
    if nnz(predictorNames(scalarRows) == scalarLabels(scalarIndex)) > 1
        scalarLabels(scalarIndex) = sprintf( ...
            "%s (%s)", ...
            scalarLabels(scalarIndex), ...
            predictorTypes(infoRow));
    end
end

%% Determine one common neuron ordering

switch options.sortMode

    case "cluster"

        neuronOrder = clusterNeuronProfiles( ...
            K, ...
            kernelRows, ...
            scalarRows);

    case "peakLag"

        if nKernelGroups == 0
            warning( ...
                "plotKernelHeatmaps:NoKernelForPeakSort", ...
                "No lagged predictor is available; retaining neuron order.");

            neuronOrder = (1:nSelectedNeurons).';

        else
            referenceKernel = kernelValues{1};

            finiteReference = referenceKernel;
            finiteReference(~isfinite(finiteReference)) = 0;

            [~,peakIndex] = max(abs(finiteReference),[],2);
            [~,neuronOrder] = sort(peakIndex,"ascend");
        end

    case "kernelStrength"

        if nKernelGroups == 0
            warning( ...
                "plotKernelHeatmaps:NoKernelForStrengthSort", ...
                "No lagged predictor is available; retaining neuron order.");

            neuronOrder = (1:nSelectedNeurons).';

        else
            referenceKernel = kernelValues{1};
            referenceKernel(~isfinite(referenceKernel)) = 0;

            kernelStrength = vecnorm(referenceKernel,2,2);
            [~,neuronOrder] = sort(kernelStrength,"descend");
        end

    case "none"

        neuronOrder = (1:nSelectedNeurons).';
end

neuronOrder = neuronOrder(:);

%% Determine shared symmetric colour limits

allDisplayedValues = [];

for kernelIndex = 1:nKernelGroups
    values = kernelValues{kernelIndex};
    allDisplayedValues = [ ...
        allDisplayedValues; ...
        values(:)]; %#ok<AGROW>
end

for scalarIndex = 1:nScalarPredictors
    values = scalarValues{scalarIndex};
    allDisplayedValues = [ ...
        allDisplayedValues; ...
        values(:)]; %#ok<AGROW>
end

finiteValues = allDisplayedValues(isfinite(allDisplayedValues));

if isempty(finiteValues)
    colorLimit = 1;
else
    colorLimit = prctile( ...
        abs(finiteValues), ...
        options.percentileLimit);
end

if isempty(colorLimit) || ~isfinite(colorLimit) || colorLimit == 0
    colorLimit = 1;
end

%% Create variable-width tiled layout

totalTileColumns = ...
    nKernelGroups*options.kernelTileWidth + ...
    nScalarPredictors*options.scalarTileWidth;

figureWidth = max( ...
    options.minimumFigureWidth, ...
    nKernelGroups*options.kernelPanelWidth + ...
    nScalarPredictors*options.scalarPanelWidth + ...
    150);

fig = figure( ...
    "Color","w", ...
    "Position",[100 100 figureWidth options.figureHeight], ...
    "Name","Predictor kernels and weights");

layout = tiledlayout( ...
    fig, ...
    1, ...
    totalTileColumns, ...
    "TileSpacing","loose", ...
    "Padding","compact");

axesHandles = gobjects(nPanels,1);
panelIndex = 0;

%% Plot lagged predictor kernels

for kernelIndex = 1:nKernelGroups

    panelIndex = panelIndex+1;

    ax = nexttile( ...
        layout, ...
        [1 options.kernelTileWidth]);

    axesHandles(panelIndex) = ax;

    kernel = kernelValues{kernelIndex};
    lags = kernelLags{kernelIndex};

    imagesc( ...
        ax, ...
        lags, ...
        1:nSelectedNeurons, ...
        kernel(neuronOrder,:));

    clim(ax,[-colorLimit,colorLimit]);
    set(ax,"YDir","reverse");

    hold(ax,"on");

    if min(lags) <= 0 && max(lags) >= 0
        xline( ...
            ax, ...
            0, ...
            "k--", ...
            "LineWidth",1, ...
            "HandleVisibility","off");
    end

    xlabel(ax,"Lag (s)");

    kernelLabel = replace(kernelNames(kernelIndex),"_"," ");
    if strlength(kernelLabel) > 9, kernelLabel = extractBefore(kernelLabel, 10); end
    title( ...
        ax, ...
        kernelLabel, ...
        "Interpreter","none");

    formatNeuronAxis( ...
        ax, ...
        panelIndex, ...
        selectedNeuronIndices, ...
        neuronOrder, ...
        options.sortMode, ...
        options.showNeuronNumbers);
end

%% Plot epoch and continuous predictor weights

for scalarIndex = 1:nScalarPredictors

    panelIndex = panelIndex+1;

    ax = nexttile( ...
        layout, ...
        [1 options.scalarTileWidth]);

    axesHandles(panelIndex) = ax;

    values = scalarValues{scalarIndex};

    imagesc( ...
        ax, ...
        1, ...
        1:nSelectedNeurons, ...
        values(neuronOrder));

    clim(ax,[-colorLimit,colorLimit]);
    set(ax,"YDir","reverse");

    xlim(ax,[0.5 1.5]);
    xticks(ax,[]);
    xlabel(ax,"");

    title( ...
        ax, ...
        replace(scalarLabels(scalarIndex),"_"," "), ...
        "Interpreter","none");

    formatNeuronAxis( ...
        ax, ...
        panelIndex, ...
        selectedNeuronIndices, ...
        neuronOrder, ...
        options.sortMode, ...
        options.showNeuronNumbers);
end

%% Figure formatting

colormap(fig,options.colormapName);

% All panels use identical neuron rows.
linkaxes(axesHandles,"y");

set( ...
    axesHandles, ...
    "YLim",[0.5 nSelectedNeurons+0.5], ...
    "YLimMode","manual");

% Shared colourbar beside the tiled layout.
cb = colorbar(axesHandles(end));

if isprop(cb,"Layout")
    cb.Layout.Tile = "east";
end

if options.standardize
    cb.Label.String = "Standardized coefficient";
else
    cb.Label.String = "Coefficient";
end

% %reduce size of colorbar
% pos = cb.Position;
% original_height = pos(4);
% new_height = original_height * 0.2;
% height_difference = original_height - new_height;
% pos(4) = new_height;
% pos(2) = pos(2) + (height_difference / 2);
% cb.Position = pos;

title( ...
    layout, ...
    sprintf( ...
        "Predictor kernels & weights: %d/%d neurons with DE >= %d%%", ...
        nnz(neuronMask), ...
        numel(neuronMask), ...
        100*options.minimumDevianceExplained), ...
    "Interpreter","none");

%% Store plotting information in the figure

setappdata( ...
    fig, ...
    "selectedNeuronIndices", ...
    selectedNeuronIndices);

setappdata( ...
    fig, ...
    "neuronOrderWithinSelection", ...
    neuronOrder);

setappdata( ...
    fig, ...
    "orderedNeuronIndices", ...
    selectedNeuronIndices(neuronOrder));

setappdata( ...
    fig, ...
    "kernelNames", ...
    kernelNames);

setappdata( ...
    fig, ...
    "scalarPredictorRows", ...
    scalarRows);

end


function neuronOrder = clusterNeuronProfiles(K,kernelRows,scalarRows)
%CLUSTERNEURONPROFILES Order neurons by similarity of predictor profiles.
%
% Each predictor group is approximately equally weighted. A lagged
% predictor with many time bins is divided by sqrt(number of bins), so it
% does not dominate simply because it has more columns.
%
% Individual coefficient dimensions are subsequently z-scored across
% neurons before clustering.

nNeurons = size(K,2);
nKernelGroups = numel(kernelRows);

featureBlocks = cell(nKernelGroups+numel(scalarRows),1);
blockIndex = 0;

%% Lagged predictor blocks

for kernelIndex = 1:nKernelGroups

    rows = kernelRows{kernelIndex};

    block = K(rows,:).';

    % Give each predictor group approximately equal overall influence,
    % regardless of the number of lag bins it contains.
    block = block ./ sqrt(max(numel(rows),1));

    blockIndex = blockIndex+1;
    featureBlocks{blockIndex} = block;
end

%% Scalar predictor blocks

for scalarIndex = 1:numel(scalarRows)

    row = scalarRows(scalarIndex);

    blockIndex = blockIndex+1;
    featureBlocks{blockIndex} = K(row,:).';
end

featureMatrix = horzcat(featureBlocks{1:blockIndex});

%% Replace nonfinite entries with the feature median

for featureIndex = 1:size(featureMatrix,2)

    values = featureMatrix(:,featureIndex);
    finiteMask = isfinite(values);

    if any(finiteMask)
        replacementValue = median(values(finiteMask));
    else
        replacementValue = 0;
    end

    values(~finiteMask) = replacementValue;
    featureMatrix(:,featureIndex) = values;
end

%% Standardize each coefficient dimension across neurons

featureMean = mean(featureMatrix,1);
featureStd = std(featureMatrix,0,1);

validFeatures = ...
    isfinite(featureStd) & ...
    featureStd > eps;

featureMatrix = featureMatrix(:,validFeatures);

if isempty(featureMatrix) || nNeurons < 2
    neuronOrder = (1:nNeurons).';
    return
end

featureMatrix = ...
    (featureMatrix-mean(featureMatrix,1)) ./ ...
    std(featureMatrix,0,1);

featureMatrix(~isfinite(featureMatrix)) = 0;

%% Hierarchical clustering, with an SVD fallback

try

    pairwiseDistance = pdist(featureMatrix,"euclidean");

    if isempty(pairwiseDistance) || ...
            all(pairwiseDistance == 0)

        neuronOrder = (1:nNeurons).';
        return
    end

    linkageTree = linkage(pairwiseDistance,"average");

    % optimalleaforder tends to place similar neighbouring clusters closer
    % together and produces cleaner-looking heatmap blocks.
    try
        neuronOrder = optimalleaforder( ...
            linkageTree, ...
            pairwiseDistance);
    catch
        [~,~,neuronOrder] = dendrogram( ...
            linkageTree, ...
            0);
        close(gcf);
    end

catch clusteringError

    warning( ...
        "plotKernelHeatmaps:ClusteringFallback", ...
        "Hierarchical clustering could not be performed (%s). " + ...
         "Falling back to ordering along the first principal axis.", ...
        clusteringError.message);

    centeredFeatures = ...
        featureMatrix-mean(featureMatrix,1);

    [~,~,rightVectors] = svd( ...
        centeredFeatures, ...
        "econ");

    firstAxisScores = ...
        centeredFeatures*rightVectors(:,1);

    [~,neuronOrder] = sort(firstAxisScores,"ascend");
end

neuronOrder = neuronOrder(:);

end


function formatNeuronAxis( ...
    ax, panelIndex, selectedNeuronIndices, neuronOrder, ...
    sortMode, showNeuronNumbers)
%FORMATNEURONAXIS Apply common neuron-axis formatting.

nNeurons = numel(neuronOrder);

ylim(ax,[0.5 nNeurons+0.5]);

sortString = "";
    switch sortMode
        case "cluster"
            sortString = " (sorted with hierarchical clustering)";
        case "peakLag"
            sortString = " (sorted by peak lag)";
        case "kernelStrength"
            sortString = " (sorted by kernel strength)";
    end

if showNeuronNumbers && ...
        nNeurons <= 50 && ...
        panelIndex == 1

    yticks(ax,1:nNeurons);

    yticklabels( ...
        ax, ...
        string(selectedNeuronIndices(neuronOrder)));

    ylabel(ax,"Neuron"+sortString);

elseif panelIndex==1

    yticks(ax,[]);
    ylabel(ax,"Neuron"+sortString);

else

    yticks(ax,[]);
    ylabel(ax,"");
end

box(ax,"off");
ax.TickDir = "out";
ax.Layer = "top";

end