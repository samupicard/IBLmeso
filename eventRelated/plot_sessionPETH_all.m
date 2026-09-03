function out = plot_sessionPETH_all(datpath, varargin)
% plot_sessionPETH_all
%
% Scans a session folder for available PETH types and plots one difference
% heatmap per type:
%   DIFF = mean(last half of conditions) - mean(first half of conditions)
%
% Preferred plotting order:
%   1) stimOn_contrastDiff
%   2) choiceMovement_choice
%   3) feedback_feedbackType
%   4) stimOn_probabilityLeft
%
% Key behavior:
%   - discovers PETH types automatically from session folder
%   - each PETH type uses its own time base
%   - types with missing/inconsistent time base are skipped with warning
%   - ROIs are filtered per FOV
%   - only ROIs present in ALL plotted PETH types are retained
%   - subsampling is performed once on that shared ROI set
%   - sortReferenceSubtype changes only the row order, not ROI membership
%   - optional onlyChronic filter using per-FOV clusterUID loading
%   - optional per-ROI normalization across all plotted types
%
% Requirements:
%   readNPY
%   brewermap
%   read_uid_csv
%   loadStructureTree / regionTokenMap if using region filtering

%% parse inputs
p = inputParser;

p.addParameter('classifierThresh', 0.5, @(x)isnumeric(x) && isscalar(x));
p.addParameter('onlyResponsive', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('onlyTuned', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('onlyChronic', false, @(x)islogical(x) || isnumeric(x));

p.addParameter('region', {}, @(x) iscell(x) && all(cellfun(@ischar, x)));
p.addParameter('stPath', 'C:\Users\Samuel\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv', @(s)ischar(s) || isstring(s));

% Can be:
%   []        -> use get_twins(evnt, subtype) for each type
%   [t0 t1]   -> same window for all types
%   struct    -> fields named with matlab.lang.makeValidName(pethType)
p.addParameter('twin_ev', [], @(x)isnumeric(x) || isstruct(x) || isempty(x));

p.addParameter('sortAbs', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('sortReferenceSubtype', '', @(x)ischar(x) || isstring(x));

p.addParameter('normalizeRows', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('caxis', [], @(x)isnumeric(x) && (isempty(x) || numel(x)==2));
p.addParameter('showColorbar', true, @(x)islogical(x) || isnumeric(x));

p.addParameter('subsampleN', [], @(x)isnumeric(x) && (isempty(x) || (isscalar(x) && x>0)));
p.addParameter('subsampleSeed', 0, @(x)isnumeric(x) && isscalar(x));

p.addParameter('pethTypes', {}, @(x) iscell(x) || isstring(x));
p.addParameter('excludePethTypes', {}, @(x) iscell(x) || isstring(x));

p.parse(varargin{:});
opt = p.Results;

opt.stPath = char(opt.stPath);
opt.sortAbs = logical(opt.sortAbs);
opt.normalizeRows = logical(opt.normalizeRows);
opt.showColorbar = logical(opt.showColorbar);
opt.onlyChronic = logical(opt.onlyChronic);
opt.sortReferenceSubtype = char(string(opt.sortReferenceSubtype));

if ~isempty(opt.subsampleN)
    opt.subsampleN = round(opt.subsampleN);
end

%% locate FOV folders
alfDir = fullfile(datpath, 'alf');
if ~exist(alfDir, 'dir')
    error('Cannot find alf folder: %s', alfDir);
end

d = dir(fullfile(alfDir, 'FOV_*'));
d = d([d.isdir]);
if isempty(d)
    error('No FOV_* folders found in %s', alfDir);
end

% deterministic FOV ordering
fovNums = nan(numel(d),1);
for i = 1:numel(d)
    tok = regexp(d(i).name, 'FOV_(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        fovNums(i) = str2double(tok{1});
    end
end
[~, ord] = sort(fovNums);
d = d(ord);

%% discover PETH types
preferredOrder = { ...
    'stimOn_contrastDiff', ...
    'choiceMovement_choice', ...
    'feedback_feedbackType', ...
    'stimOn_probabilityLeft'};

pethTypes = discoverPethTypesAcrossSession(d);

if isempty(pethTypes)
    error('No PETH types found in %s', alfDir);
end

if ~isempty(opt.pethTypes)
    req = cellstr(string(opt.pethTypes));
    pethTypes = pethTypes(ismember(pethTypes, req));
end

if ~isempty(opt.excludePethTypes)
    exc = cellstr(string(opt.excludePethTypes));
    pethTypes = pethTypes(~ismember(pethTypes, exc));
end

if isempty(pethTypes)
    error('No PETH types left after include/exclude filtering.');
end

ordered = {};
for i = 1:numel(preferredOrder)
    if ismember(preferredOrder{i}, pethTypes)
        ordered{end+1} = preferredOrder{i}; %#ok<AGROW>
    end
end
extras = pethTypes(~ismember(pethTypes, preferredOrder));
pethTypes = [ordered, extras];

%% validate time base for each type
validType = false(1, numel(pethTypes));
T_byType = cell(1, numel(pethTypes));

for s = 1:numel(pethTypes)
    pethType = pethTypes{s};
    [evnt, subtype] = splitPethType(pethType);

    ok = true;
    T_ref = [];

    for i = 1:numel(d)
        fovFolder = fullfile(d(i).folder, d(i).name);

        oddFile  = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s_odd.npy',  evnt, subtype));
        evenFile = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s_even.npy', evnt, subtype));
        avgFile  = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s.npy',      evnt, subtype));
        timeFile = fullfile(fovFolder, sprintf('PETHavgNorm_%s_%s.timeValues.npy',    evnt, subtype));

        hasPeth = (isfile(oddFile) && isfile(evenFile)) || isfile(avgFile);
        if ~hasPeth
            continue
        end

        if ~isfile(timeFile)
            warning('Skipping %s: missing timeValues file in %s', pethType, fovFolder);
            ok = false;
            break
        end

        try
            T_here = readNPY(timeFile);
            T_here = T_here(:)';
        catch ME
            warning('Skipping %s: failed to load timeValues in %s (%s)', pethType, fovFolder, ME.message);
            ok = false;
            break
        end

        if isempty(T_ref)
            T_ref = T_here;
        else
            if numel(T_ref) ~= numel(T_here) || any(abs(T_ref - T_here) > 1e-12)
                warning('Skipping %s: inconsistent time base across FOVs', pethType);
                ok = false;
                break
            end
        end
    end

    if ok && ~isempty(T_ref)
        validType(s) = true;
        T_byType{s} = T_ref;
    elseif ok && isempty(T_ref)
        warning('Skipping %s: no matching files found in any FOV', pethType);
    end
end

pethTypes = pethTypes(validType);
T_byType = T_byType(validType);

if isempty(pethTypes)
    error('No valid PETH types remained after time-base validation.');
end

nTypes = numel(pethTypes);

%% optional region setup
regionIdsExpanded = [];
if ~isempty(opt.region)
    st = loadStructureTree(opt.stPath);

    acr_all = string(st.acronym);
    rid_all = double(st.id);

    if ismember('depth', st.Properties.VariableNames)
        depth_all = double(st.depth);
    else
        error('Structure tree table must contain a "depth" column.');
    end

    regList = {};
    for r = 1:numel(opt.region)
        tok = opt.region{r};
        acrList = regionTokenMap(tok);
        if ~isempty(acrList)
            regList = [regList, acrList(:)']; %#ok<AGROW>
        else
            regList = [regList, {tok}]; %#ok<AGROW>
        end
    end

    for r = 1:numel(regList)
        reg = char(regList{r});

        if endsWith(reg, '*')
            prefix = reg(1:end-1);
            match = startsWith(acr_all, string(prefix));
            regionIdsExpanded = [regionIdsExpanded; rid_all(match)]; %#ok<AGROW>
            continue
        end

        exact = (acr_all == string(reg));
        if ~any(exact)
            continue
        end

        regionIdsExpanded = [regionIdsExpanded; rid_all(exact)]; %#ok<AGROW>

        dpth = depth_all(find(exact,1,'first'));
        if dpth >= 7
            match = startsWith(acr_all, string(reg));
            regionIdsExpanded = [regionIdsExpanded; rid_all(match)]; %#ok<AGROW>
        end
    end

    regionIdsExpanded = unique(regionIdsExpanded);

    if isempty(regionIdsExpanded)
        warning('None of the requested regions/tokens matched: %s', strjoin(opt.region, ', '));
    end

    if numel(opt.region) == 1
        regs = opt.region{1};
    else
        regs = strjoin(opt.region, ', ');
    end
else
    regs = 'All';
end

%% initialize per-type storage
Diff_even = cell(1, nTypes);
Diff_odd = cell(1, nTypes);
sortMetric = cell(1, nTypes);
globalROI_byType = cell(1, nTypes);
fov_byType = cell(1, nTypes);

for s = 1:nTypes
    Diff_even{s} = [];
    Diff_odd{s} = [];
    sortMetric{s} = [];
    globalROI_byType{s} = [];
    fov_byType{s} = [];
end

% shared metadata for kept ROIs
cellScore_all = [];
brainId_all = [];
fovLabel_all = [];
isChronic_all = [];
globalROI_all = [];

nextGlobalROI = 0;
warnedFallback = false;

%% per-FOV loading
for i = 1:numel(d)
    fovFolder = fullfile(d(i).folder, d(i).name);

    cellFile = fullfile(fovFolder, 'mpciROIs.cellClassifier.npy');
    brnFile  = fullfile(fovFolder, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy');
    typeFile = fullfile(fovFolder, 'mpciROIs.mpciROITypes.npy');
    respFile = fullfile(fovFolder, 'mpciROIs.taskResponsiveP.tsv');
    tunedFile = fullfile(fovFolder, 'mpciROIs.taskTunedP.tsv');
    cuidFile = fullfile(fovFolder, 'mpciROIs.clusterUIDs.csv');

    if ~isfile(cellFile) || ~isfile(brnFile) || ~isfile(typeFile)
        warning('Missing classifier/brainLocation/roiType files in %s (skipping FOV)', fovFolder);
        continue
    end

    cellScore = readNPY(cellFile);  cellScore = cellScore(:);
    brainIds  = readNPY(brnFile);   brainIds  = brainIds(:);
    roiTypes  = readNPY(typeFile);  roiTypes  = roiTypes(:);

    if numel(cellScore) ~= numel(brainIds) || numel(cellScore) ~= numel(roiTypes)
        warning('ROI metadata size mismatch in %s (skipping FOV)', fovFolder);
        continue
    end

    nRois = numel(cellScore);

    % somatic only
    keepType = (double(roiTypes) == 1);
    if ~any(keepType)
        continue
    end

    % chronic mask
    if opt.onlyChronic
        if isfile(cuidFile)
            try
                cUIDs = string(read_uid_csv(cuidFile));
                cUIDs = cUIDs(:);

                if numel(cUIDs) ~= nRois
                    warning('clusterUID count mismatch in %s (expected %d, got %d). Treating all ROIs as non-chronic.', ...
                        fovFolder, nRois, numel(cUIDs));
                    cUIDs = strings(nRois,1);
                end
            catch ME
                warning('Could not load %s with read_uid_csv (%s). Treating all ROIs as non-chronic.', ...
                    cuidFile, ME.message);
                cUIDs = strings(nRois,1);
            end
        else
            cUIDs = strings(nRois,1);
        end
        chronicMask = strlength(strtrim(cUIDs)) > 0;
    else
        chronicMask = true(nRois,1);
    end

    % responsive filter
    taskResp = true(nRois,1);
    if opt.onlyResponsive
        if isfile(respFile)
            try
                respTab = readtable(respFile, ...
                    'FileType','text', ...
                    'Delimiter','\t', ...
                    'VariableNamingRule','preserve');
                Resp = table2array(respTab);

                if size(Resp,1) ~= nRois
                    warning('Responsive table size mismatch in %s; using all ROIs.', fovFolder);
                else
                    taskResp = any(Resp < (0.001 / size(Resp,2)), 2);
                end
            catch ME
                warning('Could not load %s (%s). Using all ROIs for responsive filter.', respFile, ME.message);
            end
        else
            warning('Could not find mpciROIs.taskResponsiveP.tsv in %s; using all ROIs.', fovFolder);
        end
    end

    % tuned filter
    taskTuned = true(nRois,1);
    if opt.onlyTuned
        if isfile(tunedFile)
            try
                tunedTab = readtable(tunedFile, ...
                    'FileType','text', ...
                    'Delimiter','\t', ...
                    'VariableNamingRule','preserve');
                tunedNames = string(tunedTab.Properties.VariableNames);
                tunedPvals = table2array(tunedTab);

                if size(tunedPvals,1) ~= nRois
                    warning('Tuned table size mismatch in %s; using all ROIs.', fovFolder);
                else
                    if islogical(opt.onlyTuned)
                        taskTuned = any(tunedPvals < (0.05 / size(tunedPvals,2)), 2);
                    else
                        selIdx = [];

                        if isnumeric(opt.onlyTuned)
                            selIdx = opt.onlyTuned(:)';
                        elseif ischar(opt.onlyTuned) || isstring(opt.onlyTuned)
                            selNames = string(opt.onlyTuned);
                            selIdx = find(ismember(tunedNames, selNames));
                        elseif iscellstr(opt.onlyTuned) || iscell(opt.onlyTuned)
                            selNames = string(opt.onlyTuned);
                            selIdx = find(ismember(tunedNames, selNames));
                        end

                        selIdx = selIdx(selIdx >= 1 & selIdx <= size(tunedPvals,2));

                        if isempty(selIdx)
                            warning('Requested tuning selection not found in %s; using all ROIs.', fovFolder);
                        else
                            selP = tunedPvals(:, selIdx);
                            taskTuned = any((selP < 0.025) | (selP > 0.975), 2);
                        end
                    end
                end
            catch ME
                warning('Could not load %s (%s). Using all ROIs for tuning filter.', tunedFile, ME.message);
            end
        else
            warning('Could not find mpciROIs.taskTunedP.tsv in %s; using all ROIs.', fovFolder);
        end
    end

    % region filter
    if ~isempty(opt.region)
        if isempty(regionIdsExpanded)
            isRegion = false(nRois,1);
        else
            isRegion = ismember(double(brainIds), regionIdsExpanded);
        end
    else
        isRegion = true(nRois,1);
    end

    % final FOV ROI mask
    keepROI = keepType(:) & ...
              (cellScore(:) > opt.classifierThresh) & ...
              logical(taskResp(:)) & ...
              logical(taskTuned(:)) & ...
              logical(isRegion(:));

    if opt.onlyChronic
        keepROI = keepROI & logical(chronicMask(:));
    end

    if ~any(keepROI)
        continue
    end

    nKeep = sum(keepROI);
    globalROI_thisFOV = nextGlobalROI + (1:nKeep)';
    nextGlobalROI = nextGlobalROI + nKeep;

    tok = regexp(d(i).name, 'FOV_(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        fovNum = str2double(tok{1});
    else
        fovNum = nan;
    end

    % append shared kept metadata
    cellScore_all = cat(1, cellScore_all, cellScore(keepROI));
    brainId_all   = cat(1, brainId_all, brainIds(keepROI));
    fovLabel_all  = cat(1, fovLabel_all, repmat(fovNum, nKeep, 1));
    isChronic_all = cat(1, isChronic_all, chronicMask(keepROI));
    globalROI_all = cat(1, globalROI_all, globalROI_thisFOV);

    % append each type present in this FOV using same kept ROI set
    for s = 1:nTypes
        pethType = pethTypes{s};
        [evnt, subtype] = splitPethType(pethType);
        T = T_byType{s};

        oddFile  = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s_odd.npy',  evnt, subtype));
        evenFile = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s_even.npy', evnt, subtype));
        avgFile  = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s.npy',      evnt, subtype));

        hasType = false;
        Podd = [];
        Peven = [];

        if isfile(oddFile) && isfile(evenFile)
            Podd = readNPY(oddFile);
            Peven = readNPY(evenFile);
            hasType = true;
        elseif isfile(avgFile)
            Podd = readNPY(avgFile);
            Peven = Podd;
            hasType = true;
            if ~warnedFallback
                warning('Odd/Even files not found for at least one type. Falling back to avgNorm for both odd/even.');
                warnedFallback = true;
            end
        end

        if ~hasType
            continue
        end

        if size(Podd,1) ~= nRois || size(Peven,1) ~= nRois
            warning('Skipping %s in %s: ROI count mismatch.', pethType, fovFolder);
            continue
        end
        if size(Podd,3) ~= numel(T) || size(Peven,3) ~= numel(T)
            warning('Skipping %s in %s: time dimension mismatch.', pethType, fovFolder);
            continue
        end
        if size(Podd,2) < 2 || size(Peven,2) < 2
            warning('Skipping %s in %s: fewer than 2 conditions.', pethType, fovFolder);
            continue
        end

        Podd = Podd(keepROI,:,:);
        Peven = Peven(keepROI,:,:);

        nConds = size(Peven,2);
        nCondsToAvg = floor(nConds/2);
        iConds1 = 1:nCondsToAvg;
        iConds2 = nConds - fliplr(iConds1) + 1;

        diffEven_here = squeeze(mean(Peven(:,iConds2,:), 2, 'omitmissing') - ...
                                mean(Peven(:,iConds1,:), 2, 'omitmissing'));
        diffOdd_here  = squeeze(mean(Podd(:,iConds2,:), 2, 'omitmissing') - ...
                                mean(Podd(:,iConds1,:), 2, 'omitmissing'));

        % flip sign for contrastDiff only
        if strcmp(subtype, 'contrastDiff')
            diffEven_here = -diffEven_here;
            diffOdd_here  = -diffOdd_here;
        end

        twin_here = resolveTwinEv(opt.twin_ev, pethType, evnt, subtype);
        if isempty(twin_here) || numel(twin_here) ~= 2
            warning('Skipping %s in %s: invalid twin_ev.', pethType, fovFolder);
            continue
        end

        tMask = (T >= twin_here(1)) & (T <= twin_here(2));
        if ~any(tMask)
            warning('Skipping %s in %s: twin_ev [%g %g] does not overlap time base.', ...
                pethType, fovFolder, twin_here(1), twin_here(2));
            continue
        end

        dResp_here = mean(diffOdd_here(:,tMask), 2);

        Diff_even{s} = cat(1, Diff_even{s}, diffEven_here);
        Diff_odd{s}  = cat(1, Diff_odd{s}, diffOdd_here);
        sortMetric{s} = cat(1, sortMetric{s}, dResp_here);
        globalROI_byType{s} = cat(1, globalROI_byType{s}, globalROI_thisFOV);
        fov_byType{s} = cat(1, fov_byType{s}, repmat(fovNum, nKeep, 1));
    end
end

if isempty(cellScore_all)
    error('No ROIs remained after metadata filtering.');
end

%% keep usable types only
hasRows = cellfun(@(x) ~isempty(x), Diff_even);
pethTypes = pethTypes(hasRows);
T_byType = T_byType(hasRows);
Diff_even = Diff_even(hasRows);
Diff_odd = Diff_odd(hasRows);
sortMetric = sortMetric(hasRows);
globalROI_byType = globalROI_byType(hasRows);
fov_byType = fov_byType(hasRows);
nTypes = numel(pethTypes);

if nTypes == 0
    error('No PETH types remained after loading and twin_ev checks.');
end

%% restrict to ROIs shared across all plotted types
sharedIDs = globalROI_byType{1};
for s = 2:nTypes
    sharedIDs = intersect(sharedIDs, globalROI_byType{s}, 'stable');
end

if isempty(sharedIDs)
    error('No ROIs are shared across all plotted PETH types after filtering.');
end

for s = 1:nTypes
    localMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for k = 1:numel(globalROI_byType{s})
        localMap(double(globalROI_byType{s}(k))) = double(k);
    end

    localRows = zeros(numel(sharedIDs),1);
    for k = 1:numel(sharedIDs)
        localRows(k) = localMap(double(sharedIDs(k)));
    end

    Diff_even{s} = Diff_even{s}(localRows,:);
    Diff_odd{s} = Diff_odd{s}(localRows,:);
    sortMetric{s} = sortMetric{s}(localRows);
    globalROI_byType{s} = globalROI_byType{s}(localRows);
    fov_byType{s} = fov_byType{s}(localRows);
end

% reduce shared metadata to shared IDs too
metaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
for k = 1:numel(globalROI_all)
    metaMap(double(globalROI_all(k))) = double(k);
end
metaRows = zeros(numel(sharedIDs),1);
for k = 1:numel(sharedIDs)
    metaRows(k) = metaMap(double(sharedIDs(k)));
end

cellScore_all = cellScore_all(metaRows);
brainId_all = brainId_all(metaRows);
fovLabel_all = fovLabel_all(metaRows);
isChronic_all = isChronic_all(metaRows);
globalROI_all = globalROI_all(metaRows);

%% optional row normalization across all plotted types using shared global ROI ids
roiScale = [];
if opt.normalizeRows
    roiScaleMap = containers.Map('KeyType', 'double', 'ValueType', 'double');

    for s = 1:nTypes
        ids = globalROI_byType{s};
        vals = max(abs(Diff_even{s}), [], 2);

        for k = 1:numel(ids)
            id = double(ids(k));
            v = double(vals(k));

            if isKey(roiScaleMap, id)
                roiScaleMap(id) = max(roiScaleMap(id), v);
            else
                roiScaleMap(id) = v;
            end
        end
    end

    roiScale = zeros(numel(sharedIDs),1);
    for k = 1:numel(sharedIDs)
        roiScale(k) = roiScaleMap(double(sharedIDs(k)));
    end
    roiScale(~isfinite(roiScale) | roiScale == 0) = 1;

    for s = 1:nTypes
        Diff_even{s} = Diff_even{s} ./ roiScale;
    end
end

%% shared subsampling on common ROI set
if ~isempty(opt.subsampleN)
    nAvail = numel(sharedIDs);

    if opt.subsampleN < nAvail
        sel = chooseSubsample(nAvail, opt.subsampleN, opt.subsampleSeed);

        sharedIDs = sharedIDs(sel);
        if ~isempty(roiScale)
            roiScale = roiScale(sel);
        end

        for s = 1:nTypes
            Diff_even{s} = Diff_even{s}(sel,:);
            Diff_odd{s} = Diff_odd{s}(sel,:);
            sortMetric{s} = sortMetric{s}(sel);
            globalROI_byType{s} = globalROI_byType{s}(sel);
            fov_byType{s} = fov_byType{s}(sel);
        end

        cellScore_all = cellScore_all(sel);
        brainId_all = brainId_all(sel);
        fovLabel_all = fovLabel_all(sel);
        isChronic_all = isChronic_all(sel);
        globalROI_all = globalROI_all(sel);
    end
end

%% sorting
sharedOrderMode = false;
sharedOrderType = '';
sortIdxByType = cell(1, nTypes);

if ~isempty(opt.sortReferenceSubtype)
    ref = char(string(opt.sortReferenceSubtype));
    refIdx = find(strcmp(pethTypes, ref), 1);

    if isempty(refIdx)
        warning('sortReferenceSubtype "%s" not found among plotted types. Using per-type sorting.', ref);

        for s = 1:nTypes
            if opt.sortAbs
                [~, sortIdxByType{s}] = sort(abs(sortMetric{s}), 'descend');
            else
                [~, sortIdxByType{s}] = sort(sortMetric{s}, 'descend');
            end
        end
    else
        if opt.sortAbs
            [~, refOrder] = sort(abs(sortMetric{refIdx}), 'descend');
        else
            [~, refOrder] = sort(sortMetric{refIdx}, 'descend');
        end

        for s = 1:nTypes
            sortIdxByType{s} = refOrder;
        end

        sharedOrderMode = true;
        sharedOrderType = ref;
    end
else
    for s = 1:nTypes
        if opt.sortAbs
            [~, sortIdxByType{s}] = sort(abs(sortMetric{s}), 'descend');
        else
            [~, sortIdxByType{s}] = sort(sortMetric{s}, 'descend');
        end
    end
end

%% prepare images for plotting
Diff_plot = cell(1, nTypes);
for s = 1:nTypes
    Diff_plot{s} = Diff_even{s}(sortIdxByType{s}, :);
end

%% caxis
if isempty(opt.caxis)
    if opt.normalizeRows
        opt.caxis = [-1 1];
    else
        allVals = [];
        for s = 1:nTypes
            allVals = [allVals; Diff_plot{s}(:)]; %#ok<AGROW>
        end

        clim = prctile(allVals, [1 99.9]);
        lim = round(max(abs(clim)),2);
        if ~isfinite(lim) || lim == 0
            lim = max(abs(allVals));
        end
        if ~isfinite(lim) || lim == 0
            lim = 1;
        end
        opt.caxis = [-lim lim];
    end
end

%% plot
figW = max(200, 20 + 80*nTypes);
figH = 500;

figure('Position', [1500 -100 figW figH], ...
       'Color', 'w', ...
       'Name', sprintf('%s | all session PETHs', sessNameFromPath(datpath)));

tl = tiledlayout(1, nTypes, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

tl.Units = 'normalized';
tl.OuterPosition = [0.02 0.02 0.85 0.85];

ax = gobjects(1, nTypes);

for s = 1:nTypes
    ax(s) = nexttile;
    T = T_byType{s};

    imagesc(T, 1:size(Diff_plot{s},1), Diff_plot{s});
    colormap(ax(s), brewermap([], '*RdBu'));
    caxis(ax(s), opt.caxis);
    axis(ax(s), 'tight');

    set(ax(s), ...
        'XTick', [], ...
        'YTick', [], ...
        'XColor', 'none', ...
        'YColor', 'none', ...
        'TickLength', [0 0], ...
        'Box', 'off');

    ttl = prettyPethTitle(pethTypes{s});
    title(ax(s), ttl, 'Interpreter', 'none');

    hold(ax(s), 'on');
    xline(ax(s), 0, 'k-');
    hold(ax(s), 'off');
end

if opt.showColorbar
    drawnow;

    refPos = ax(end).Position;
    cbH = refPos(4) / 10;
    cbW = 0.02;
    cbX = 0.83;
    cbY = refPos(2) + (refPos(4) - cbH)/2;

    axCB = axes('Position', [cbX cbY cbW cbH], ...
                'Visible', 'off', ...
                'Color', 'none');

    colormap(axCB, brewermap([], '*RdBu'));
    caxis(axCB, opt.caxis);

    cb = colorbar(axCB, 'Position', [cbX cbY cbW cbH]);
    cb.YAxisLocation = 'right';
    cb.TickDirection = 'out';
    cb.Box = 'off';
    cb.FontSize = 9;
    cb.Ticks = linspace(opt.caxis(1), opt.caxis(2), 3);

    if opt.normalizeRows
        cb.Label.String = 'Normalized \Delta response';
    else
        cb.Label.String = '\Delta response';
    end
    cb.Label.Interpreter = 'tex';
    cb.Label.Rotation = 270;
    cb.Label.VerticalAlignment = 'bottom';
end

drawnow;
axBL = ax(1);
pos = axBL.Position;

dx_sec = 1;
nPlot = size(Diff_plot{1},1);
dy_nrn = 20 * max(1, round(max(1,nPlot)/500));
fprintf('scalebar: x = %.1f seconds, y = %d ROIs\n', dx_sec, dy_nrn);

T0 = T_byType{1};
Tspan = T0(end) - T0(1);
if Tspan <= 0
    Tspan = 1;
end
if nPlot <= 0
    nPlot = 1;
end

dx_norm = (dx_sec / Tspan) * pos(3);
dy_norm = (dy_nrn / nPlot) * pos(4);

marginX = 0.010;
marginY = 0.025;
x0 = max(pos(1) - marginX, 0.02);
y0 = max(pos(2) - marginY, 0.02);

annotation(gcf, 'line', [x0, x0 + dx_norm], [y0, y0], ...
    'Color', 'k', 'LineWidth', 2);
annotation(gcf, 'line', [x0, x0], [y0, y0 + dy_norm], ...
    'Color', 'k', 'LineWidth', 2);

% ----- title text -----
nShown = size(Diff_plot{1}, 1);

topLine = sprintf('%s (n=%d)', regs, nShown);
midLine = '\Delta resp. for even trials';

if ~isempty(opt.sortReferenceSubtype)
    ttl = prettyPethTitle(opt.sortReferenceSubtype);
    sortLbl = ttl{1}; % reuse top line (e.g. 'STIM SIDE')
    botLine = sprintf('sorted by odd trial \\Delta resp. for %s', sortLbl);
else
    botLine = 'sorted by odd \Delta resp.';
end

sgtitle(sprintf('%s\n%s\n%s', topLine, midLine, botLine), ...
    'Interpreter', 'tex', ...
    'FontSize', 10, ...
    'FontWeight', 'bold');

% move sgtitle upward
drawnow;
sg = findall(gcf, 'Type', 'Text', 'Tag', 'sgtitle');
if ~isempty(sg)
    sg.Units = 'normalized';
    sg.Position(2) = 0.975;
end

%% outputs
out = struct();
out.datpath = datpath;
out.pethTypes = pethTypes;
out.T_byType = T_byType;

out.cellScore = cellScore_all;
out.brainIds = brainId_all;
out.fovLabel = fovLabel_all;
out.isChronic = isChronic_all;

out.region = opt.region;
out.regionIdsExpanded = regionIdsExpanded;

out.Diff_even = Diff_even;
out.Diff_odd = Diff_odd;
out.Diff_plot = Diff_plot;
out.sortMetric = sortMetric;
out.sortIdxByType = sortIdxByType;

out.normalizeRows = opt.normalizeRows;
out.roiScale = roiScale;
out.sharedROI = sharedIDs;
out.globalROI_byType = globalROI_byType;

out.sortReferenceSubtype = opt.sortReferenceSubtype;
out.sharedOrderMode = sharedOrderMode;
out.sharedOrderType = sharedOrderType;

end

%% ---------------- helper functions ----------------

function pethTypes = discoverPethTypesAcrossSession(d)
pethTypes = {};

for i = 1:numel(d)
    fovFolder = fullfile(d(i).folder, d(i).name);
    files = dir(fullfile(fovFolder, 'mpciROIs.PETHavgNorm_*.npy'));

    for k = 1:numel(files)
        nm = files(k).name;

        tok = regexp(nm, '^mpciROIs\.PETHavgNorm_(.+?)(?:_(odd|even))?\.npy$', 'tokens', 'once');
        if isempty(tok)
            continue
        end

        core = tok{1};
        us = strfind(core, '_');
        if isempty(us)
            continue
        end

        pethTypes{end+1} = core; %#ok<AGROW>
    end
end

pethTypes = unique(pethTypes, 'stable');
end

function [evnt, subtype] = splitPethType(pethType)
parts = split(string(pethType), "_");
if numel(parts) < 2
    error('Invalid pethType: %s', pethType);
end
evnt = char(parts(1));
subtype = char(strjoin(parts(2:end), "_"));
end

function twin = resolveTwinEv(twinOpt, pethType, evnt, subtype)
if isempty(twinOpt)
    [defaultTwinEv, ~, ~] = get_twins(evnt, subtype);
    twin = defaultTwinEv;
    return
end

if isnumeric(twinOpt)
    twin = twinOpt;
    return
end

if isstruct(twinOpt)
    fn = matlab.lang.makeValidName(pethType);
    if isfield(twinOpt, fn)
        twin = twinOpt.(fn);
    else
        [defaultTwinEv, ~, ~] = get_twins(evnt, subtype);
        twin = defaultTwinEv;
    end
    return
end

twin = [];
end

function sel = chooseSubsample(nAvail, nTake, seed)
if seed ~= 0
    rng(seed);
    sel = sort(randperm(nAvail, nTake));
else
    sel = unique(round(linspace(1, nAvail, nTake)));
    while numel(sel) < nTake
        cand = setdiff(1:nAvail, sel);
        sel = sort([sel, cand(1)]);
    end
end
end

function nm = sessNameFromPath(datpath)
sp = split(string(datpath), filesep);
if numel(sp) >= 3
    nm = sprintf('%s | %s | %s', sp(end-2), sp(end-1), sp(end));
else
    nm = char(datpath);
end
end

function ttl = prettyPethTitle(pethType)
% Return a 2-line title as {topLine; bottomLine}

switch char(string(pethType))
    case 'stimOn_contrastDiff'
        ttl = {'STIM SIDE'; 'stimOn'};
    case 'choiceMovement_choice'
        ttl = {'CHOICE'; 'choiceMov'};
    case 'feedback_feedbackType'
        ttl = {'FEEDBACK'; 'feedback'};
    case 'stimOn_probabilityLeft'
        ttl = {'BLOCK'; 'stimOn'};
    otherwise
        % fallback: subtype on top, event on bottom
        [evnt, subtype] = splitPethType(pethType);

        % optional light prettifying for fallback
        topLine = upper(regexprep(subtype, '([a-z])([A-Z])', '$1 $2'));
        botLine = evnt;

        ttl = {topLine; botLine};
    end
end