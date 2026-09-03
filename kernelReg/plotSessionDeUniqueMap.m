function plotSessionDeUniqueMap( ...
    sessionPaths, predictorNames, options)
%PLOTSESSIONDEUNIQUEMAP Plot deviance explained across one or more sessions.
%
% Data from every FOV in every supplied session are concatenated into the
% same anatomical maps.
%
% One subplot is produced for each predictor. Each subplot shows:
%
%       x     = mediolateral coordinate
%       y     = anteroposterior coordinate
%       color = unique cross-validated deviance explained
%
% An additional subplot shows full-model deviance explained.
%
% Only ROIs with mpciROITypes == options.roiTypeValue are included.
%
% Expected files within each FOV folder:
%
%   mpciROIs.mlapdv_estimate.npy
%   mpciROIs.mpciROITypes.npy
%   mpciROIs.deFull.npy
%   mpciROIs.deUnique_<predictor>.npy
%
% INPUTS
% ------
% sessionPaths
%   One session path or multiple session paths, supplied as a string
%   vector, character vector, or cell array of character vectors.
%
% predictorNames
%   Optional predictor names corresponding to the saved deUnique files.
%   If empty, predictor names are inferred across all discovered FOVs.
%
% OPTIONS
% -------
% markerSize
%   Scatter marker area. Default: 5.
%
% markerAlpha
%   Scatter marker opacity. Default: 0.5.
%
% nColors
%   Number of colormap entries. Default: 256.
%
% colorLimit
%   Symmetric unique-deviance color limit in percent. When empty, it is
%   determined from all finite unique-deviance values.
%
% sharedColorLimits
%   Use the same unique-deviance color limits for every predictor.
%
% minimumDE
%   Optional minimum full-model deviance explained, in percent.
%
% fullColorLimit
%   Upper color limit for full-model deviance explained, in percent.
%
% roiTypeValue
%   ROI type interpreted as a cell. Default: 1.
%
% figurePosition
%   Figure position.
%
% bas
%   Top-down atlas structure returned by aratopdown.atlas.build_topdown.
%
% REQUIREMENTS
% ------------
%   readNPY
%   brewermap

arguments
    sessionPaths
    predictorNames (1,:) string = strings(1,0)

    options.markerSize (1,1) double ...
        {mustBePositive} = 4

    options.markerAlpha (1,1) double ...
        {mustBeGreaterThanOrEqual(options.markerAlpha,0), ...
         mustBeLessThanOrEqual(options.markerAlpha,1)} = 0.4

    options.nColors (1,1) double ...
        {mustBeInteger,mustBePositive} = 256

    options.colorLimit double = 10
    options.sharedColorLimits (1,1) logical = true

    options.minimumDE double = 2

    options.fullColorLimit (1,1) double ...
        {mustBePositive} = 20

    options.roiTypeValue (1,1) double = 1

    options.figurePosition (1,4) double = ...
        [100 50 900 800]

    options.bas (1,1) struct = ...
        aratopdown.atlas.build_topdown
end

%% Normalize inputs

sessionPaths = normalizeSessionPaths(sessionPaths);
predictorNames = predictorNames(:).';

if isempty(sessionPaths)
    error( ...
        "plotSessionDeUniqueMap:NoSessions", ...
        "At least one session path must be supplied.");
end

if ~isempty(options.minimumDE) && ...
        (~isscalar(options.minimumDE) || ...
         ~isfinite(options.minimumDE))

    error( ...
        "plotSessionDeUniqueMap:InvalidMinimumDE", ...
        "options.minimumDE must be empty or a finite scalar in percent.");
end

%% Validate session folders and discover FOVs

nRequestedSessions = numel(sessionPaths);

allRoiTypeFiles = struct( ...
    "name",{}, ...
    "folder",{}, ...
    "date",{}, ...
    "bytes",{}, ...
    "isdir",{}, ...
    "datenum",{}, ...
    "sessionIndex",{});

sessionHasFOVs = false(nRequestedSessions,1);

for sessionIndex = 1:nRequestedSessions

    thisSessionPath = sessionPaths(sessionIndex);

    if ~isfolder(thisSessionPath)
        warning( ...
            "plotSessionDeUniqueMap:SessionNotFound", ...
            "Skipping missing session folder:\n%s", ...
            thisSessionPath);
        continue
    end

    thisRoiTypeFiles = dir(fullfile( ...
        thisSessionPath, ...
        "**", ...
        "mpciROIs.mpciROITypes.npy"));

    if isempty(thisRoiTypeFiles)
        warning( ...
            "plotSessionDeUniqueMap:NoFOVsInSession", ...
            "No FOV folders were found beneath:\n%s", ...
            thisSessionPath);
        continue
    end

    sessionHasFOVs(sessionIndex) = true;

    for fileIndex = 1:numel(thisRoiTypeFiles)

        thisFile = thisRoiTypeFiles(fileIndex);
        thisFile.sessionIndex = sessionIndex;

        allRoiTypeFiles(end+1) = thisFile; %#ok<AGROW>
    end
end

includedSessionPaths = sessionPaths(sessionHasFOVs);

if isempty(allRoiTypeFiles)
    error( ...
        "plotSessionDeUniqueMap:NoFOVs", ...
        "No FOV folders containing 'mpciROIs.mpciROITypes.npy' were found in any supplied session.");
end

nSessions = numel(includedSessionPaths);
nFOVs = numel(allRoiTypeFiles);

%% Infer predictor names when not supplied

if isempty(predictorNames)

    discoveredNames = strings(0,1);

    for fovIndex = 1:nFOVs

        fovPath = string(allRoiTypeFiles(fovIndex).folder);

        deUniqueFiles = dir(fullfile( ...
            fovPath, ...
            "mpciROIs.deUnique_*.npy"));

        for fileIndex = 1:numel(deUniqueFiles)

            fileName = string(deUniqueFiles(fileIndex).name);

            predictorName = erase( ...
                fileName, ...
                "mpciROIs.deUnique_");

            predictorName = erase( ...
                predictorName, ...
                ".npy");

            discoveredNames(end+1,1) = ...
                predictorName; %#ok<AGROW>
        end
    end

    predictorNames = unique( ...
        discoveredNames, ...
        "stable").';

    if isempty(predictorNames)
        error( ...
            "plotSessionDeUniqueMap:NoPredictors", ...
            "No files matching 'mpciROIs.deUnique_*.npy' were found in any supplied session.");
    end
end

%% Apply preferred predictor order

preferredOrder = [ ...
    "goCue", ...
    "stimulusLeft", ...
    "stimulusRight", ...
    "choice", ...
    "action", ...
    "feedback", ...
    "block", ...
    "wheelSpeed"];

orderedPreferred = preferredOrder( ...
    ismember(preferredOrder,predictorNames));

remainingPredictors = predictorNames( ...
    ~ismember(predictorNames,preferredOrder));

predictorNames = [ ...
    orderedPreferred, ...
    remainingPredictors];

nPredictors = numel(predictorNames);

fprintf( ...
    "Aggregating %d sessions, %d FOV folders and %d predictors.\n", ...
    nSessions, ...
    nFOVs, ...
    nPredictors);

%% Load and concatenate FOV data

allML = cell(nFOVs,1);
allAP = cell(nFOVs,1);
allDeFull = cell(nFOVs,1);
allSessionIndex = cell(nFOVs,1);

deUniqueByPredictor = cell(nPredictors,1);

for predictorIndex = 1:nPredictors
    deUniqueByPredictor{predictorIndex} = cell(nFOVs,1);
end

for fovIndex = 1:nFOVs

    fovPath = string(allRoiTypeFiles(fovIndex).folder);

    sourceSessionIndex = ...
        allRoiTypeFiles(fovIndex).sessionIndex;

    sourceSessionPath = ...
        sessionPaths(sourceSessionIndex);

    fprintf( ...
        "Loading FOV %d/%d: %s\n", ...
        fovIndex, ...
        nFOVs, ...
        fovPath);

    coordinateFile = fullfile( ...
        fovPath, ...
        "mpciROIs.mlapdv_estimate.npy");

    roiTypeFile = fullfile( ...
        fovPath, ...
        "mpciROIs.mpciROITypes.npy");

    deFullFile = fullfile( ...
        fovPath, ...
        "mpciROIs.deFull.npy");

    if ~isfile(deFullFile)
        error( ...
            "plotSessionDeUniqueMap:MissingDeFull", ...
            ["The full-model deviance file is missing from:\n%s\n" ...
             "Session: %s"], ...
            fovPath, ...
            sourceSessionPath);
    end

    if ~isfile(coordinateFile)
        error( ...
            "plotSessionDeUniqueMap:MissingCoordinates", ...
            ["Coordinate file is missing from FOV folder:\n%s\n" ...
             "Session: %s"], ...
            fovPath, ...
            sourceSessionPath);
    end

    coordinates = double(readNPY(coordinateFile));

    roiTypes = double(readNPY(roiTypeFile));
    roiTypes = roiTypes(:);

    deFull = double(readNPY(deFullFile));
    deFull = deFull(:);

    nROIs = numel(roiTypes);

    if size(coordinates,2) < 2
        error( ...
            "plotSessionDeUniqueMap:InvalidCoordinates", ...
            ["Coordinates in '%s' must contain at least two columns " ...
             "for ML and AP."], ...
            coordinateFile);
    end

    if size(coordinates,1) ~= nROIs
        error( ...
            "plotSessionDeUniqueMap:ROICountMismatch", ...
            ["Coordinate and ROI-type files have different numbers " ...
             "of ROIs in:\n%s"], ...
            fovPath);
    end

    if numel(deFull) ~= nROIs
        error( ...
            "plotSessionDeUniqueMap:DeFullCountMismatch", ...
            ["The deFull file contains %d values, but the FOV " ...
             "contains %d ROIs:\n%s"], ...
            numel(deFull), ...
            nROIs, ...
            fovPath);
    end

    isCell = ...
        roiTypes == options.roiTypeValue;

    validCoordinate = ...
        isfinite(coordinates(:,1)) & ...
        isfinite(coordinates(:,2));

    validDeFull = isfinite(deFull);

    % deFull is stored as a fraction; filtering is specified in percent.
    if isempty(options.minimumDE)
        passesMinimumDE = true(size(deFull));
    else
        passesMinimumDE = ...
            100*deFull >= options.minimumDE;
    end

    keepROI = ...
        isCell & ...
        validCoordinate & ...
        validDeFull & ...
        passesMinimumDE;

    allML{fovIndex} = ...
        coordinates(keepROI,1);

    allAP{fovIndex} = ...
        coordinates(keepROI,2);

    allDeFull{fovIndex} = ...
        deFull(keepROI);

    allSessionIndex{fovIndex} = repmat( ...
        sourceSessionIndex, ...
        nnz(keepROI), ...
        1);

    for predictorIndex = 1:nPredictors

        predictorName = ...
            predictorNames(predictorIndex);

        deUniqueFile = fullfile( ...
            fovPath, ...
            sprintf( ...
                "mpciROIs.deUnique_%s.npy", ...
                predictorName));

        if ~isfile(deUniqueFile)
            error( ...
                "plotSessionDeUniqueMap:MissingDeUnique", ...
                ["The deUnique file for predictor '%s' is missing " ...
                 "from:\n%s\nSession: %s"], ...
                predictorName, ...
                fovPath, ...
                sourceSessionPath);
        end

        deUnique = double(readNPY(deUniqueFile));
        deUnique = deUnique(:);

        if numel(deUnique) ~= nROIs
            error( ...
                "plotSessionDeUniqueMap:DeUniqueCountMismatch", ...
                ["The deUnique file for predictor '%s' contains %d " ...
                 "values, but the FOV contains %d ROIs:\n%s"], ...
                predictorName, ...
                numel(deUnique), ...
                nROIs, ...
                fovPath);
        end

        deUniqueByPredictor{predictorIndex}{fovIndex} = ...
            deUnique(keepROI);
    end
end

%% Concatenate across FOVs and sessions

ml = vertcat(allML{:});
ap = vertcat(allAP{:});

deFull = vertcat(allDeFull{:});
deFullPercent = 100*deFull;

cellSessionIndex = vertcat(allSessionIndex{:});

for predictorIndex = 1:nPredictors

    deUniqueByPredictor{predictorIndex} = ...
        vertcat( ...
            deUniqueByPredictor{predictorIndex}{:});
end

if isempty(ml)
    error( ...
        "plotSessionDeUniqueMap:NoCellROIs", ...
        ["No finite ROIs with mpciROITypes == %g passed the " ...
         "requested filtering across the supplied sessions."], ...
        options.roiTypeValue);
end

if isempty(options.minimumDE)

    fprintf( ...
        "Plotting %d cells across %d sessions and %d FOVs.\n", ...
        numel(ml), ...
        nSessions, ...
        nFOVs);

else

    fprintf( ...
        "Plotting %d cells across %d sessions and %d FOVs with " +...
         "full-model deviance explained >= %.2f%%.\n", ...
        numel(ml), ...
        nSessions, ...
        nFOVs, ...
        options.minimumDE);
end

%% Determine color limits

if ~isempty(options.colorLimit)

    if ~isscalar(options.colorLimit) || ...
            ~isfinite(options.colorLimit) || ...
            options.colorLimit <= 0

        error( ...
            "plotSessionDeUniqueMap:InvalidColorLimit", ...
            "options.colorLimit must be empty or a positive scalar.");
    end

    sharedColorLimit = ...
        options.colorLimit;

else

    allFiniteValues = [];

    for predictorIndex = 1:nPredictors

        % Convert to percent before determining limits.
        values = ...
            100*deUniqueByPredictor{predictorIndex};

        allFiniteValues = [ ...
            allFiniteValues; ...
            values(isfinite(values))]; %#ok<AGROW>
    end

    if isempty(allFiniteValues)
        sharedColorLimit = 1;
    else
        sharedColorLimit = ...
            max(abs(allFiniteValues));
    end

    if sharedColorLimit == 0
        sharedColorLimit = 1;
    end
end

%% Set common spatial limits

mlLimits = paddedLimits(ml);
apLimits = paddedLimits(ap);

%% Create figure

nPanels = nPredictors; %nPanels = nPredictors+1;

%nTileColumns = 3; nTileRows = ceil(nPanels/nTileColumns);
nTileColumns = 8; nTileRows = 1;

figureWidth = max( ...
    options.figurePosition(3), ...
    300*nTileColumns);

figureHeight = max( ...
    options.figurePosition(4), ...
    260*nTileRows);

fig = figure( ...
    "Color","w", ...
    "Position",[ ...
        options.figurePosition(1), ...
        options.figurePosition(2), ...
        figureWidth, ...
        figureHeight], ...
    "Name","Aggregated deUnique maps");

layout = tiledlayout( ...
    fig, ...
    nTileRows, ...
    nTileColumns, ...
    "TileSpacing","tight", ...
    "Padding","tight");

uniqueColorMap = ...
    brewermap(options.nColors,'*RdBu');

baseColorMap = ...
    brewermap(256,'*RdBu');

redHalf = ...
    baseColorMap(129:256,:);

redColorMap = interp1( ...
    linspace(0,1,size(redHalf,1)), ...
    redHalf, ...
    linspace(0,1,options.nColors), ...
    "linear");

axesHandles = gobjects(1,nPanels);

%% Plot predictor-specific unique deviance

for predictorIndex = 1:nPredictors

    ax = nexttile(layout,predictorIndex);
    axesHandles(predictorIndex) = ax;

    hold(ax,"on");

    deUniquePercent = ...
        100*deUniqueByPredictor{predictorIndex};

    validValue = ...
        isfinite(deUniquePercent);

    scatter( ...
        ax, ...
        ml(validValue), ...
        ap(validValue), ...
        options.markerSize, ...
        deUniquePercent(validValue), ...
        "filled", ...
        "MarkerFaceAlpha",options.markerAlpha, ...
        "MarkerEdgeAlpha",options.markerAlpha);

    if options.sharedColorLimits

        thisColorLimit = ...
            sharedColorLimit;

    else

        finiteValues = ...
            deUniquePercent(isfinite(deUniquePercent));

        if isempty(finiteValues)
            thisColorLimit = 1;
        else
            thisColorLimit = ...
                max(abs(finiteValues));
        end

        if thisColorLimit == 0
            thisColorLimit = 1;
        end
    end

    colormap(ax,uniqueColorMap);
    clim(ax,[-thisColorLimit,thisColorLimit]);

    plotAtlasBoundaries( ...
        ax, ...
        options.bas);

    formatSpatialAxes( ...
        ax, ...
        mlLimits, ...
        apLimits);

    title( ...
        ax, ...
        replace( ...
            predictorNames(predictorIndex), ...
            "_", ...
            " "), ...
        "Interpreter","none","FontSize",12);

    % Use a single unique-deviance colorbar when limits are shared.
    if ~options.sharedColorLimits || ...
            predictorIndex == nPredictors

        colorBar = colorbar(ax);

        colorBar.Label.String = ...
            "Unique c.v. deviance explained (%)";

        colorBar.Label.FontSize = 10;
        colorBar.FontSize = 10;

        colorBar.Limits = ...
            [-thisColorLimit,thisColorLimit];
    end
end

%% Plot full-model deviance explained

% axFull = nexttile(layout,nPanels);
% axesHandles(nPanels) = axFull;
% 
% hold(axFull,"on");
% 
% validFull = ...
%     isfinite(deFullPercent);
% 
% scatter( ...
%     axFull, ...
%     ml(validFull), ...
%     ap(validFull), ...
%     options.markerSize, ...
%     deFullPercent(validFull), ...
%     "filled", ...
%     "MarkerFaceAlpha",options.markerAlpha, ...
%     "MarkerEdgeAlpha",options.markerAlpha);
% 
% colormap(axFull,redColorMap);
% clim(axFull,[0,options.fullColorLimit]);
% 
% plotAtlasBoundaries( ...
%     axFull, ...
%     options.bas);
% 
% formatSpatialAxes( ...
%     axFull, ...
%     mlLimits, ...
%     apLimits);
% 
% title( ...
%     axFull, ...
%     "Full model", ...
%     "Interpreter","none");
% 
% fullColorBar = colorbar(axFull);
% 
% fullColorBar.Label.String = ...
%     "Full c.v. deviance explained (%)";

%% Link spatial axes

linkaxes(axesHandles,"xy");

set( ...
    axesHandles, ...
    "XLim",mlLimits, ...
    "YLim",apLimits, ...
    "XLimMode","manual", ...
    "YLimMode","manual");

%% Figure title

if nSessions == 1

    sourceText = includedSessionPaths(1);

else

    animals = strings(nSessions,1);

    for sessionIndex = 1:nSessions
        animals(sessionIndex) = ...
            getAnimalIdFromSessionPath( ...
                includedSessionPaths(sessionIndex));
    end

    nAnimals = numel(unique(animals));

    sourceText = sprintf( ...
        "%d sessions from %d animals", ...
        nSessions, ...
        nAnimals);
end

title( ...
    layout, ...
    sprintf( ...
        "%s\nExplainable deviance by ROI (n = %d cells, %d FOVs)", ...
        sourceText, ...
        numel(ml), ...
        nFOVs), ...
    "Interpreter","none");

%% Store aggregated source information

setappdata( ...
    fig, ...
    "sessionPaths", ...
    includedSessionPaths);

setappdata( ...
    fig, ...
    "cellSessionIndex", ...
    cellSessionIndex);

setappdata( ...
    fig, ...
    "predictorNames", ...
    predictorNames);

end


function sessionPaths = normalizeSessionPaths(inputPaths)
%NORMALIZESESSIONPATHS Convert supported path inputs to a string column.

if ischar(inputPaths)

    sessionPaths = string(inputPaths);

elseif isstring(inputPaths)

    sessionPaths = inputPaths;

elseif iscell(inputPaths)

    try
        sessionPaths = string(inputPaths);
    catch
        error( ...
            "plotSessionDeUniqueMap:InvalidSessionPaths", ...
            ["sessionPaths must be a path, string vector, or cell " ...
             "array containing paths."]);
    end

else

    error( ...
        "plotSessionDeUniqueMap:InvalidSessionPaths", ...
        ["sessionPaths must be a character vector, string vector, " ...
         "or cell array containing paths."]);
end

sessionPaths = sessionPaths(:);
sessionPaths = strip(sessionPaths);
sessionPaths = sessionPaths(sessionPaths ~= "");

% Avoid loading the same session twice.
sessionPaths = unique(sessionPaths,"stable");

end


function animalId = getAnimalIdFromSessionPath(sessionPath)
%GETANIMALIDFROMSESSIONPATH Extract animal from animal/date/session path.

parts = regexp( ...
    char(sessionPath), ...
    "[\\/]+", ...
    "split");

parts = parts(~cellfun("isempty",parts));

if numel(parts) >= 3
    animalId = string(parts{end-2});
else
    animalId = "unknown";
end

end


function limits = paddedLimits(values)
%PADDEDLIMITS Return plotting limits with a small spatial margin.

minimumValue = min(values,[],"omitnan");
maximumValue = max(values,[],"omitnan");

valueRange = maximumValue-minimumValue;

if valueRange == 0
    padding = max(abs(minimumValue)*0.05,1);
else
    padding = 0.03*valueRange;
end

limits = [ ...
    minimumValue-padding, ...
    maximumValue+padding];

end


function plotAtlasBoundaries(ax,bas)
%PLOTATLASBOUNDARIES Add dorsal atlas boundaries to one axes.

if isempty(bas) || ...
        ~isfield(bas,"dorsal_brain_areas")
    return
end

hold(ax,"on");

brainAreas = bas.dorsal_brain_areas;

% Preserve the exclusion used in the original function, but avoid invalid
% indexing when fewer than 12 regions are present.
lastArea = max(numel(brainAreas)-11,0);

for areaIndex = 1:lastArea

    boundaries = ...
        brainAreas(areaIndex).boundaries_stereotax;

    for boundaryIndex = 1:numel(boundaries)

        boundaryCoordinates = ...
            boundaries{boundaryIndex};

        plot( ...
            ax, ...
            1000*boundaryCoordinates(:,2), ...
            1000*boundaryCoordinates(:,1), ...
            "Color",[0.6 0.6 0.6], ...
            "HandleVisibility","off");
    end
end

end


function formatSpatialAxes(ax,mlLimits,apLimits)
%FORMATSPATIALAXES Apply common anatomical map formatting.

xlim(ax,mlLimits);
ylim(ax,apLimits);

daspect(ax,[1 1 1]);

% Reinstate exact limits after setting the data aspect ratio.
ax.XLim = mlLimits;
ax.YLim = apLimits;
ax.XLimMode = "manual";
ax.YLimMode = "manual";

box(ax,"off");
axis(ax,"off");

end