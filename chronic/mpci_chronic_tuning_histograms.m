function out = mpci_chronic_tuning_histograms(sessionPaths, opt)
% MPCI chronic tuning histograms:
%   1) Loads taskTunedStat.npy (stats in [-1, 1]) instead of taskTunedP.tsv
%   3) Extracts full distributions of stats per day × region × test
%   4) For each test, plots a dates-by-region grid of HISTOGRAMS
%      - shared x/y axes and shared bin edges
%      - vertical line at 0 and at the histogram mean

if nargin < 2, opt = struct(); end

% optional region filter inputs
if ~isfield(opt,'region'), opt.region =  {'VIS*','RSP*','SS*','AUD*','MO*'}; end % acronyms/tokens
if ~isfield(opt,'st'),     opt.st = ''; end
if ~isfield(opt,'stPath'), opt.stPath = "C:\Users\Samuel\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv"; end

% histogram options (shared defaults)
if ~isfield(opt,'nBins') || isempty(opt.nBins), opt.nBins = 20; end
if ~isfield(opt,'statRange') || isempty(opt.statRange), opt.statRange = [-1 1]; end
if ~isfield(opt,'maxCols') || isempty(opt.maxCols), opt.maxCols = 18; end  % optional: cap cols shown at once

sessionPaths = string(sessionPaths(:));
if isempty(sessionPaths)
    error('sessionPaths is empty.');
end

records = table();
allTestNames = string.empty(1,0);   % union of all test columns encountered

for i = 1:numel(sessionPaths)
    sp = sessionPaths(i);

    % Accept either session root or alf path
    if endsWith(sp, filesep + "alf") || endsWith(sp, "/alf") || endsWith(sp, "\alf")
        alfPath = sp;
        sessRoot = fileparts(alfPath);
    else
        sessRoot = sp;
        alfPath = fullfile(sessRoot, 'alf');
    end
    if ~isfolder(alfPath), continue; end

    % Derive day + sessionID from path .../<day>/<sessionID>[/alf]
    [parentDir, sessionID] = fileparts(sessRoot);
    [subjDir, dayStr] = fileparts(parentDir);
    [~,subjStr] = fileparts(subjDir);
    dayStr = string(dayStr);
    sessionID = string(sessionID);

    % Enumerate FOVs
    fovs = dir(fullfile(alfPath, 'FOV*'));
    fovs = fovs([fovs.isdir]);

    for f = 1:numel(fovs)
        fovName = string(fovs(f).name);
        fovPath = fullfile(alfPath, fovName);

        % --- UPDATED: stats file(s) ---
        % Try a few common filenames for robustness
        statFileA = fullfile(fovPath, 'mpciROIs.taskTunedStat.npy');
        statFileB = fullfile(fovPath, 'taskTunedStat.npy');

        % Optional: if you also have a TSV with names, we can use it
        statTsvA  = fullfile(fovPath, 'mpciROIs.taskTunedStat.tsv');
        statTsvB  = fullfile(fovPath, 'taskTunedStat.tsv');

        regFile = fullfile(fovPath, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy');
        uidFile = fullfile(fovPath, 'mpciROIs.clusterUIDs.csv');

        % Always require stats + region
        hasStat = isfile(statFileA) || isfile(statFileB) || isfile(statTsvA) || isfile(statTsvB);
        if ~(hasStat && isfile(regFile))
            continue
        end

        % -------- Load stats + test names --------
        S = [];
        testNames = string.empty(1,0);

        if isfile(statTsvA) || isfile(statTsvB)
            % If a TSV exists, prefer it (has variable names)
            if isfile(statTsvA), tsvPath = statTsvA; else, tsvPath = statTsvB; end
            sTab = readtable(tsvPath, 'FileType','text', 'Delimiter','\t', 'VariableNamingRule','preserve');
            testNames = string(sTab.Properties.VariableNames);
            S = table2array(sTab);
        else
            % Otherwise load NPY (values in [-1,1])
            if isfile(statFileA), npyPath = statFileA; else, npyPath = statFileB; end
            S = double(readNPY(npyPath));
            if ndims(S) ~= 2
                error('Expected a 2D array in %s (nROIs x nTests).', npyPath);
            end
            % Create generic test names if none available
            testNames = "Test" + string(1:size(S,2));
        end

        n = size(S,1);

        % Update union of all tests; if new tests appear, add NaN columns to existing records
        newTests = setdiff(testNames, allTestNames, 'stable');
        if ~isempty(newTests)
            for nt = 1:numel(newTests)
                records.(newTests(nt)) = nan(height(records), 1);
            end
            allTestNames = [allTestNames, newTests];
        end

        % Load regions
        regionIds = double(readNPY(regFile));
        regionIds = regionIds(:);

        if numel(regionIds) ~= n
            warning('Size mismatch (region vs stats) in %s; skipping FOV.', fovPath);
            continue
        end

        % Load UIDs if present; otherwise fill empty strings
        if isfile(uidFile)
            uidCol = read_uid_csv(uidFile);
            if numel(uidCol) ~= n
                warning('Size mismatch (UID vs stats) in %s; filling empty UIDs for this FOV.', fovPath);
                uidCol = repmat("", n, 1);
            end
        else
            uidCol = repmat("", n, 1);
        end

        % Assemble per-ROI rows
        T = table();
        T.day        = repmat(dayStr, n, 1);
        T.sessionID  = repmat(sessionID, n, 1);
        T.fov        = repmat(fovName, n, 1);
        T.clusterUID = uidCol;
        T.regionId   = regionIds;

        % Add the test columns that exist in this file
        for k = 1:numel(testNames)
            T.(testNames(k)) = S(:,k);
        end

        % ---- Robust append to records with variable columns ----
        if isempty(records)
            records = T;
            allTestNames = testNames;
        else
            allTestNames = unique([allTestNames, testNames], 'stable');

            % 1) Add columns present in T but missing in records
            missingInRecords = setdiff(T.Properties.VariableNames, records.Properties.VariableNames, 'stable');
            for m = 1:numel(missingInRecords)
                v = missingInRecords{m};
                records.(v) = nan(height(records), 1);
            end

            % 2) Add columns present in records but missing in T
            missingInT = setdiff(records.Properties.VariableNames, T.Properties.VariableNames, 'stable');
            for m = 1:numel(missingInT)
                v = missingInT{m};
                if isnumeric(records.(v))
                    T.(v) = nan(height(T), 1);
                elseif isdatetime(records.(v))
                    T.(v) = NaT(height(T), 1);
                else
                    T.(v) = repmat("", height(T), 1);
                end
            end

            % 3) Reorder T columns to match records, then concatenate
            T = T(:, records.Properties.VariableNames);
            records = [records; T];
        end
    end
end

if isempty(records)
    error('No usable FOVs found. (Need taskTunedStat.npy/tsv + brainLocation .npy at minimum)');
end

% -----------------------------
% Day ordering (YYYY-MM-DD)
% -----------------------------
records.day = datetime(records.day, 'InputFormat','yyyy-MM-dd');
days = unique(records.day);
days = sort(days);
nDays = numel(days);

% Tests
metaCols = ["day","sessionID","fov","clusterUID","regionId"];
testCols = allTestNames;
nTests = numel(testCols);

% Only rows with non-empty UID (i.e., tracked ROIs)
hasUID = records.clusterUID ~= "" & ~ismissing(records.clusterUID);
Ruid = records(hasUID, :);

if isempty(Ruid)
    error('No tracked ROIs found (no non-empty clusterUIDs.csv entries in provided sessions).');
end

% -----------------------------
% Sanity: duplicates of same UID within the same day
% -----------------------------
[Gdup, ~, ~] = findgroups(Ruid.day, Ruid.clusterUID);
cnt = splitapply(@numel, Ruid.clusterUID, Gdup);
if any(cnt > 1)
    warning('Found UID×day groups with >1 occurrence. Your chronic definition will exclude these UIDs.');
end

%% 1) Get chronic UIDs tracked on all days
cUIDs_d1 = get_clusterUIDs(sessionPaths(1));
is_tracked = get_trackedROIs(cUIDs_d1, sessionPaths(2:end));
chronicUIDs = cUIDs_d1(is_tracked);

Rc = Ruid(ismember(Ruid.clusterUID, chronicUIDs), :);

if isempty(chronicUIDs)
    warning('No chronic UIDs found under the "exactly once per day" criterion.');
end

%% 2) Stable region per chronic UID = modal region across days
chUIDs = chronicUIDs(:);
nCh = numel(chUIDs);
stableRegion = nan(nCh,1);
tieFlag = false(nCh,1);

for j = 1:nCh
    uid = chUIDs(j);
    r = Rc.regionId(Rc.clusterUID == uid);
    [stableRegion(j), tieFlag(j)] = modal_region(r);
end

uid2reg = table(chUIDs, stableRegion, tieFlag, ...
    'VariableNames', {'clusterUID','stableRegionId','tieFlag'});

if nCh > 0
    [tf, loc] = ismember(Rc.clusterUID, uid2reg.clusterUID);
    if ~all(tf), error('UID lookup failed (unexpected).'); end
    Rc.stableRegionId = uid2reg.stableRegionId(loc);

    % Optional region filter (by stable modal region per UID)
    [regionIdsExpanded, regsLabel] = expand_region_ids(opt);

    if ~isempty(regionIdsExpanded)
        keep = ismember(double(Rc.stableRegionId), double(regionIdsExpanded));
        Rc = Rc(keep, :);

        chronicUIDs = unique(Rc.clusterUID);
        uid2reg = uid2reg(ismember(uid2reg.clusterUID, chronicUIDs), :);

        if isempty(Rc)
            warning('Region filter matched no chronic UIDs. Requested: %s', regsLabel);
        end
    else
        regsLabel = "All";
    end
else
    Rc.stableRegionId = [];
    regsLabel = "All";
end

% enforce minimum N chronic UIDs per region
minChronicPerRegion = 100;
nUIDsPerReg = [];

if ~isempty(Rc)
    [Greg, regVals] = findgroups(double(Rc.stableRegionId));
    nUIDsPerReg = splitapply(@(u) numel(unique(u)), Rc.clusterUID, Greg);

    keepReg = nUIDsPerReg >= minChronicPerRegion;
    keepRegionIds = regVals(keepReg);

    Rc = Rc(ismember(double(Rc.stableRegionId), keepRegionIds), :);

    chronicUIDs = unique(Rc.clusterUID);
    uid2reg = uid2reg(ismember(uid2reg.clusterUID, chronicUIDs), :);
    nUIDsPerReg = nUIDsPerReg(keepReg);

    if isempty(keepRegionIds)
        warning('No regions have >= %d chronic UIDs after filtering.', minChronicPerRegion);
    end
end

% Regions based on stable IDs (after filtering)
if ~isempty(Rc)
    regionList = sort(unique(Rc.stableRegionId));
else
    regionList = [];
end
nRegions = numel(regionList);

%% 3) FULL DISTRIBUTIONS per day × stable region × test
% Store vectors in a cell array:
%   distStats{di, ri, ti} = [stat_1, stat_2, ...]  (NaNs removed)
distStats = cell(nDays, nRegions, nTests);
nDen = zeros(nDays, nRegions, nTests);

for di = 1:nDays
    if nCh == 0, continue; end
    Rd = Rc(Rc.day == days(di), :);

    for ri = 1:nRegions
        rID = regionList(ri);
        regMask = (Rd.stableRegionId == rID);

        if ~any(regMask), continue; end

        for ti = 1:nTests
            x = Rd.(testCols(ti));
            x = x(regMask);

            ok = ~isnan(x);
            nDen(di,ri,ti) = sum(ok);

            if any(ok)
                distStats{di,ri,ti} = x(ok);
            else
                distStats{di,ri,ti} = [];
            end
        end
    end
end

%% 4) Plot dates-by-region grid of histograms for each test

% Get structure tree for label mapping
if isfield(opt,'st') && ~isempty(opt.st)
    st = opt.st;
elseif isfield(opt,'stPath') && strlength(string(opt.stPath)) > 0
    st = loadStructureTree(opt.stPath);
else
    error('Structure tree required for region acronym labeling.');
end

rid_all = double(st.id);
acr_all = string(st.acronym);

% Map region IDs -> acronyms
regionLabels = strings(size(regionList));
region_ix = nan(size(regionList));

for ii = 1:numel(regionList)
    idx = find(rid_all == regionList(ii), 1);
    region_ix(ii) = idx;
    if ~isempty(idx)
        regionLabels(ii) = acr_all(idx-1); % hack to remove '1' (kept from your original)
    else
        regionLabels(ii) = string(regionList(ii));
    end
end

% sorting of region (by structure tree row index)
if nRegions > 0 && ~isempty(regionList)
    [~, idxRow] = ismember(double(regionList), rid_all);
    [~, sortIdx] = sort(idxRow);

    regionList   = regionList(sortIdx);
    regionLabels = regionLabels(sortIdx);

    distStats = distStats(:, sortIdx, :);
    nDen      = nDen(:, sortIdx, :);

    if ~isempty(nUIDsPerReg)
        nUIDsPerReg = nUIDsPerReg(sortIdx);
    end
end

% shared bin edges for all plots (per test)
edges = linspace(opt.statRange(1), opt.statRange(2), opt.nBins+1);

for ti = 1:nTests

    % Which days have any data for this test?
    dayHasData = squeeze(sum(nDen(:,:,ti), 2)) > 0;  % nDays x 1  (sum across regions)
    dayIdxAll = find(dayHasData);

    if isempty(dayIdxAll)
        warning('No data for test %s; skipping plot.', testCols(ti));
        continue
    end

    % Optional: show in chunks if many days
    dayChunks = chunk_indices(dayIdxAll, opt.maxCols);

    for ci = 1:numel(dayChunks)
        dayIdx = dayChunks{ci};
        nCols = numel(dayIdx);

        % Precompute max count for shared y-limits (across all dates of each region)
        maxCounts = nan(1,nRegions);
        allMeans = nan(nRegions, nCols);

        for ri = 1:nRegions
            maxCount = 0;
            for c = 1:nCols
                di = dayIdx(c);
                x = distStats{di,ri,ti};
                if isempty(x), continue; end
                counts = histcounts(x, edges);
                maxCount = max(maxCount, max(counts));
                allMeans(ri,c) = mean(x);
            end
            if maxCount == 0, maxCount = 1; end
            maxCounts(ri) = maxCount;
        end

        figName = sprintf('%s | %s | %s', subjStr, testCols(ti), regsLabel);
        figure('Color','w', 'Name', figName, 'Position', [1100, 50, 80*nCols + 200, 80*nRegions + 100]);

        tl = tiledlayout(nRegions, nCols, 'TileSpacing','compact', 'Padding','compact');

        for ri = 1:nRegions
            for c = 1:nCols
                di = dayIdx(c);
                ax = nexttile(tl, (ri-1)*nCols + c);

                x = distStats{di,ri,ti};

                if isempty(x)
                    histogram(ax, [], edges); % empty
                else
                    histogram(ax, x, edges);
                end

                xlim(ax, opt.statRange);
                ylim(ax, [0 maxCounts(ri)]);

                % Vertical lines: 0 and mean
                xline(ax, 0, 'k-', 'LineWidth', 1);
                if ~isnan(allMeans(ri,c))
                    xline(ax, allMeans(ri,c), 'r-', 'LineWidth', 1);
                end

                % Labels: top row = day, first col = region
                if ri == 1
                    title(ax, string(days(di), 'yyyy-MM-dd'), 'Interpreter','none', 'FontSize', 9);
                end
                if c == 1
                    if ~isempty(nUIDsPerReg)
                        ylabel(ax, regionLabels(ri) + " (" + string(nUIDsPerReg(ri)) + ")");
                    else
                        ylabel(ax, regionLabels(ri));
                    end
                end

                % reduce clutter
                if ri ~= nRegions
                    ax.XTickLabel = {};
                end
                if c ~= 1
                    ax.YTickLabel = {};
                end
                ax.Box = 'on';
            end
        end

        % overall labels
        xlabel(tl, 'Task-tuned statistic ([-1, 1])');
        ylabel(tl, 'Region');
        title(tl, sprintf('%s\n%s (hist grid) | days %d-%d of %d', ...
            subjStr, testCols(ti), find(dayHasData,1,'first'), find(dayHasData,1,'last'), nDays), ...
            'Interpreter','none');
    end
end

% Output
out = struct();
out.recordsAll      = records;
out.recordsTracked  = Ruid;
out.recordsChronic  = Rc;
out.days            = days;
out.regionList      = regionList;
out.regionLabels    = regionLabels;
out.testCols        = testCols;
out.distStats       = distStats;   % cell(nDays,nRegions,nTests) of vectors
out.nDen            = nDen;
out.chronicUIDs     = chronicUIDs;
out.uid2reg         = uid2reg;
out.histEdges       = edges;
out.opt             = opt;

end

%% helper: modal region
function [m, tied] = modal_region(r)
r = r(:);
r = r(~isnan(r));
if isempty(r)
    m = NaN; tied = false; return
end
u = unique(r);
counts = zeros(numel(u),1);
for k = 1:numel(u)
    counts(k) = sum(r == u(k));
end
maxc = max(counts);
modes = u(counts == maxc);
tied = numel(modes) > 1;
m = min(modes); % deterministic tie-break
end

%% helper: chunk day indices into groups of at most maxCols
function chunks = chunk_indices(idx, maxCols)
if isempty(idx)
    chunks = {};
    return
end
if isempty(maxCols) || maxCols <= 0 || numel(idx) <= maxCols
    chunks = {idx};
    return
end
n = numel(idx);
nChunks = ceil(n / maxCols);
chunks = cell(1, nChunks);
for i = 1:nChunks
    a = (i-1)*maxCols + 1;
    b = min(i*maxCols, n);
    chunks{i} = idx(a:b);
end
end