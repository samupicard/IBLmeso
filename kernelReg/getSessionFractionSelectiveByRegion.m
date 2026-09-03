function results = getSessionFractionSelectiveByRegion( ...
    sessionPath, predictorNames, options)
%GETSESSIONFRACTIONSELECTIVEBYREGION Fraction selective by Allen region.
%
% All FOVs belonging to one session are combined. Valid cell ROIs are
% assigned to Allen atlas regions using:
%
%   mpciROIs.brainLocationIds_ccf_2017_estimate.npy
%
% For each region and predictor, the fraction of selective cells is:
%
%                        number of selective cells
%   fractionSelective = ---------------------------
%                         number of valid cells
%
% A cell is selective for predictor p when:
%
%   deFull   > options.fullThresholdPercent
%   deUnique > options.uniqueThresholdPercent
%
% By default, these thresholds are:
%
%   deFull   > 2%
%   deUnique > 1%
%
% deFull and deUnique are stored as fractions in the NPY files and are
% converted to percentages internally.
%
% INPUTS
% ------
% sessionPath
%   Path to one session folder.
%
% predictorNames
%   Optional predictor names. If empty, names are inferred from files
%   matching:
%
%       mpciROIs.deUnique_*.npy
%
% OPTIONS
% -------
% minimumDE
%   Optional minimum full-model deviance explained, in percent, used to
%   filter the population before calculating regional fractions.
%
%   For example, minimumDE=2 restricts the denominator to cells with at
%   least 2% full-model deviance explained.
%
%   Empty means no additional deFull filtering.
%
% fullThresholdPercent
%   Full-model selectivity threshold in percent. Default: 2.
%
% uniqueThresholdPercent
%   Unique-deviance selectivity threshold in percent. Default: 1.
%
% roiTypeValue
%   ROI type interpreted as a cell. Default: 1.
%
% validRegionIds
%   Optional Allen region IDs to retain. Empty retains all positive,
%   finite region IDs.
%
% OUTPUT
% ------
% results.regionIds
%   R-by-1 Allen atlas region IDs.
%
% results.predictorNames
%   1-by-P predictor names.
%
% results.fractionSelective
%   R-by-P fraction selective.
%
% results.nSelective
%   R-by-P number of selective cells.
%
% results.nValid
%   R-by-P number of cells entering each fraction.
%
% results.nCellsByRegion
%   R-by-1 number of cells passing the base ROI filters.
%
% results.table
%   Long-form table with one row per region and predictor.
%
% REQUIREMENT
% -----------
%   readNPY

arguments
    sessionPath (1,1) string
    predictorNames (1,:) string = strings(1,0)

    options.minimumDE double = []
    options.fullThresholdPercent (1,1) double = 2
    options.uniqueThresholdPercent (1,1) double = 1
    options.roiTypeValue (1,1) double = 1
    options.validRegionIds (:,1) double = zeros(0,1)

    options.uniqueFractionThresholdPercent (1,1) double = 1
end

sessionPath = string(sessionPath);
predictorNames = predictorNames(:).';

%% Validate options

if ~isfolder(sessionPath)
    error( ...
        "getSessionFractionSelectiveByRegion:SessionNotFound", ...
        "Session folder does not exist: %s", ...
        sessionPath);
end

if ~isempty(options.minimumDE) && ...
        (~isscalar(options.minimumDE) || ...
         ~isfinite(options.minimumDE))

    error( ...
        "getSessionFractionSelectiveByRegion:InvalidMinimumDE", ...
        "options.minimumDE must be empty or a finite scalar in percent.");
end

%% Find FOV folders

roiTypeFiles = dir(fullfile( ...
    sessionPath, ...
    "**", ...
    "mpciROIs.mpciROITypes.npy"));

if isempty(roiTypeFiles)
    error( ...
        "getSessionFractionSelectiveByRegion:NoFOVs", ...
        "No FOV folders containing " + ...
        "'mpciROIs.mpciROITypes.npy' were found beneath:\n%s", ...
        sessionPath);
end

nFOVs = numel(roiTypeFiles);

%% Infer predictor names

if isempty(predictorNames)

    discoveredNames = strings(0,1);

    for fovIndex = 1:nFOVs

        fovPath = string(roiTypeFiles(fovIndex).folder);

        predictorFiles = dir(fullfile( ...
            fovPath, ...
            "mpciROIs.deUnique_*.npy"));

        for fileIndex = 1:numel(predictorFiles)

            fileName = string(predictorFiles(fileIndex).name);

            predictorName = erase( ...
                fileName, ...
                "mpciROIs.deUnique_");

            predictorName = erase( ...
                predictorName, ...
                ".npy");

            discoveredNames(end+1,1) = predictorName; %#ok<AGROW>
        end
    end

    predictorNames = unique(discoveredNames,"stable").';

    if isempty(predictorNames)
        error( ...
            "getSessionFractionSelectiveByRegion:NoPredictors", ...
            "No files matching 'mpciROIs.deUnique_*.npy' were " + ...
            "found beneath:\n%s", ...
            sessionPath);
    end
end

%% Apply preferred predictor order

preferredOrder = [ ...
    "stimulusLeft", ...
    "stimulusRight", ...
    "choice", ...
    "feedback", ...
    "block", ...
    "wheelSpeed"];

orderedPreferred = preferredOrder( ...
    ismember(preferredOrder,predictorNames));

remainingPredictors = predictorNames( ...
    ~ismember(predictorNames,preferredOrder));

predictorNames = [orderedPreferred,remainingPredictors];

nPredictors = numel(predictorNames);

%% Allocate FOV-level storage

regionByFOV = cell(nFOVs,1);
deFullByFOV = cell(nFOVs,1);

deUniqueByPredictor = cell(nPredictors,1);

for predictorIndex = 1:nPredictors
    deUniqueByPredictor{predictorIndex} = cell(nFOVs,1);
end

%% Load FOVs

for fovIndex = 1:nFOVs

    fovPath = string(roiTypeFiles(fovIndex).folder);

    roiTypeFile = fullfile( ...
        fovPath, ...
        "mpciROIs.mpciROITypes.npy");

    regionFile = fullfile( ...
        fovPath, ...
        "mpciROIs.brainLocationIds_ccf_2017_estimate.npy");

    deFullFile = fullfile( ...
        fovPath, ...
        "mpciROIs.deFull.npy");

    if ~isfile(regionFile)
        error( ...
            "getSessionFractionSelectiveByRegion:MissingRegionIds", ...
            "The Allen region ID file is missing from:\n%s", ...
            fovPath);
    end

    if ~isfile(deFullFile)
        error( ...
            "getSessionFractionSelectiveByRegion:MissingDeFull", ...
            "The full-model deviance file is missing from:\n%s", ...
            fovPath);
    end

    roiTypes = double(readNPY(roiTypeFile));
    regionIds = double(readNPY(regionFile));
    deFull = double(readNPY(deFullFile));

    roiTypes = roiTypes(:);
    regionIds = regionIds(:);
    deFull = deFull(:);

    nROIs = numel(roiTypes);

    if numel(regionIds) ~= nROIs
        error( ...
            "getSessionFractionSelectiveByRegion:RegionCountMismatch", ...
            "The region-ID file contains %d values, but the FOV " + ...
            "contains %d ROIs:\n%s", ...
            numel(regionIds), ...
            nROIs, ...
            fovPath);
    end

    if numel(deFull) ~= nROIs
        error( ...
            "getSessionFractionSelectiveByRegion:DeFullCountMismatch", ...
            "The deFull file contains %d values, but the FOV " + ...
            "contains %d ROIs:\n%s", ...
            numel(deFull), ...
            nROIs, ...
            fovPath);
    end

    %% Base ROI filter

    isCell = roiTypes == options.roiTypeValue;

    validRegion = ...
        isfinite(regionIds) & ...
        regionIds > 0 & ...
        regionIds == round(regionIds);

    validDeFull = isfinite(deFull);

    if isempty(options.validRegionIds)
        requestedRegion = true(nROIs,1);
    else
        requestedRegion = ismember( ...
            regionIds, ...
            options.validRegionIds);
    end

    if isempty(options.minimumDE)
        passesMinimumDE = true(nROIs,1);
    else
        passesMinimumDE = ...
            100*deFull >= options.minimumDE;
    end

    keepROI = ...
        isCell & ...
        validRegion & ...
        validDeFull & ...
        requestedRegion & ...
        passesMinimumDE;

    regionByFOV{fovIndex} = regionIds(keepROI);
    deFullByFOV{fovIndex} = deFull(keepROI);

    %% Load predictor-specific unique deviance

    for predictorIndex = 1:nPredictors

        predictorName = predictorNames(predictorIndex);

        deUniqueFile = fullfile( ...
            fovPath, ...
            sprintf( ...
                "mpciROIs.deUnique_%s.npy", ...
                predictorName));

        if ~isfile(deUniqueFile)
            error( ...
                "getSessionFractionSelectiveByRegion:MissingDeUnique", ...
                "The deUnique file for predictor '%s' is missing " + ...
                "from:\n%s", ...
                predictorName, ...
                fovPath);
        end

        deUnique = double(readNPY(deUniqueFile));
        deUnique = deUnique(:);

        if numel(deUnique) ~= nROIs
            error( ...
                "getSessionFractionSelectiveByRegion:" + ...
                "DeUniqueCountMismatch", ...
                "The deUnique file for predictor '%s' contains %d " + ...
                "values, but the FOV contains %d ROIs:\n%s", ...
                predictorName, ...
                numel(deUnique), ...
                nROIs, ...
                fovPath);
        end

        deUniqueByPredictor{predictorIndex}{fovIndex} = ...
            deUnique(keepROI);
    end
end

%% Concatenate all FOVs

regionIdsPerCell = vertcat(regionByFOV{:});
deFullPerCell = vertcat(deFullByFOV{:});

for predictorIndex = 1:nPredictors
    deUniqueByPredictor{predictorIndex} = ...
        vertcat(deUniqueByPredictor{predictorIndex}{:});
end

if isempty(regionIdsPerCell)
    error( ...
        "getSessionFractionSelectiveByRegion:NoValidCells", ...
        "No cells passed the requested session-level filters.");
end

uniqueRegionIds = unique(regionIdsPerCell);
nRegions = numel(uniqueRegionIds);

%% Calculate fractions

fractionSelective = nan(nRegions,nPredictors);
nSelective = zeros(nRegions,nPredictors);
nValid = zeros(nRegions,nPredictors);
nCellsByRegion = zeros(nRegions,1);

medianDeUnique = nan(nRegions,nPredictors);
q25DeUnique = nan(nRegions,nPredictors);
q75DeUnique = nan(nRegions,nPredictors);

fractionAboveUniqueThreshold = nan(nRegions,nPredictors);
nAboveUniqueThreshold = zeros(nRegions,nPredictors);

deFullPercent = 100*deFullPerCell;

for regionIndex = 1:nRegions

    inRegion = ...
        regionIdsPerCell == uniqueRegionIds(regionIndex);

    nCellsByRegion(regionIndex) = nnz(inRegion);

    for predictorIndex = 1:nPredictors

        deUnique = ...
            deUniqueByPredictor{predictorIndex};

        deUniquePercent = 100*deUnique;

        % A cell enters this predictor's denominator only when its
        % predictor-specific deUnique value is finite.
        validForPredictor = ...
            inRegion & ...
            isfinite(deUniquePercent);

        regionValues = deUniquePercent(validForPredictor);

        if ~isempty(regionValues)

            medianDeUnique(regionIndex,predictorIndex) = ...
                median(regionValues,"omitnan");

            q25DeUnique(regionIndex,predictorIndex) = ...
                prctile(regionValues,25);

            q75DeUnique(regionIndex,predictorIndex) = ...
                prctile(regionValues,75);

            nAboveUniqueThreshold(regionIndex,predictorIndex) = ...
                nnz(regionValues > options.uniqueThresholdPercent);

            fractionAboveUniqueThreshold(regionIndex,predictorIndex) = ...
                nAboveUniqueThreshold(regionIndex,predictorIndex) ./ ...
                nnz(regionValues > options.uniqueFractionThresholdPercent);
            
        end

        isSelective = ...
            validForPredictor & ...
            deFullPercent > options.fullThresholdPercent & ...
            deUniquePercent > options.uniqueThresholdPercent;

        nValid(regionIndex,predictorIndex) = ...
            nnz(validForPredictor);

        nSelective(regionIndex,predictorIndex) = ...
            nnz(isSelective);

        if nValid(regionIndex,predictorIndex) > 0
            fractionSelective(regionIndex,predictorIndex) = ...
                nSelective(regionIndex,predictorIndex) ./ ...
                nValid(regionIndex,predictorIndex);
        end
    end
end

%% Construct long-form output table

[regionGrid,predictorGrid] = ndgrid( ...
    uniqueRegionIds, ...
    1:nPredictors);
regionColumn = regionGrid(:);

predictorColumn = predictorNames(predictorGrid(:)).';
predictorColumn = predictorColumn(:);

resultsTable = table( ...
    regionColumn, ...
    predictorColumn, ...
    fractionSelective(:), ...
    nSelective(:), ...
    nValid(:),...
    medianDeUnique(:), ...
    q25DeUnique(:), ...
    q75DeUnique(:), ...
    fractionAboveUniqueThreshold(:));

resultsTable.Properties.VariableNames = { ...
    'regionId', ...
    'predictorName', ...
    'fractionSelective', ...
    'nSelective', ...
    'nValid', ...
    'medianDeUniquePercent', ...
    'q25DeUniquePercent', ...
    'q75DeUniquePercent', ...
    'fractionAboveUniqueThreshold'};

%% Return results

results = struct;

results.sessionPath = sessionPath;
results.predictorNames = predictorNames;
results.regionIds = uniqueRegionIds;

results.fractionSelective = fractionSelective;
results.nSelective = nSelective;
results.nValid = nValid;
results.nCellsByRegion = nCellsByRegion;

results.minimumDE = options.minimumDE;
results.fullThresholdPercent = options.fullThresholdPercent;
results.uniqueThresholdPercent = ...
    options.uniqueThresholdPercent;

results.medianDeUnique = medianDeUnique;
results.q25DeUnique = q25DeUnique;
results.q75DeUnique = q75DeUnique;

results.fractionAboveUniqueThreshold = ...
    fractionAboveUniqueThreshold;

results.nAboveUniqueThreshold = ...
    nAboveUniqueThreshold;

results.nFOVs = nFOVs;
results.nCells = numel(regionIdsPerCell);

results.table = resultsTable;

end