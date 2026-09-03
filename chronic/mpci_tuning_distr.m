function out = mpci_tuning_distr(sessionPaths, opt)
% MPCI tuning distribution summaries:
%   1) Loads taskTunedStat.tsv/npy
%   2) Extracts full distributions of stats per day × region × test
%   3) Plots compact per-region distribution summaries
%
% trackedOnly behavior:
%   opt.trackedOnly = true  -> original chronic/tracked behavior
%   opt.trackedOnly = false -> include all ROIs with mpciROITypes==1,
%                              regardless of clusterUID

if nargin < 2, opt = struct(); end

% optional region filter inputs
if ~isfield(opt,'region'), opt.region =  {'VIS*','RSP*','SS*','AUD*','MO*'}; end % acronyms/tokens
if ~isfield(opt,'st'),     opt.st = ''; end
if ~isfield(opt,'stPath'), opt.stPath = "C:\Users\Samuel\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv"; end

% plotting / stat options
if ~isfield(opt,'nBins') || isempty(opt.nBins), opt.nBins = 20; end
if ~isfield(opt,'statRange') || isempty(opt.statRange), opt.statRange = []; end
if ~isfield(opt,'violinWidth') || isempty(opt.violinWidth), opt.violinWidth = .38; end
if ~isfield(opt,'stat') || isempty(opt.stat), opt.stat = 'ccMI'; end
if ~isfield(opt,'trackedOnly') || isempty(opt.trackedOnly), opt.trackedOnly = true; end
if ~isfield(opt,'maxRegionsPerFigure') || isempty(opt.maxRegionsPerFigure), opt.maxRegionsPerFigure = 10; end
if ~isfield(opt,'minROIsPerRegion') || isempty(opt.minROIsPerRegion), opt.minROIsPerRegion = 100; end
if ~isfield(opt,'responsiveOnly') || isempty(opt.responsiveOnly), opt.responsiveOnly = false; end
if ~isfield(opt,'responsiveAlpha') || isempty(opt.responsiveAlpha), opt.responsiveAlpha = 0.005; end

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

        % stats file(s)
        statFileA = fullfile(fovPath, 'mpciROIs.taskTunedStat.npy');
        statFileB = fullfile(fovPath, 'taskTunedStat.npy');
        statTsvA  = fullfile(fovPath, 'mpciROIs.taskTunedStat.tsv');
        statTsvB  = fullfile(fovPath, 'taskTunedStat.tsv');

        regFile  = fullfile(fovPath, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy');
        uidFile  = fullfile(fovPath, 'mpciROIs.clusterUIDs.csv');
        typeFile = fullfile(fovPath, 'mpciROIs.mpciROITypes.npy');

        respFile = fullfile(fovPath, 'mpciROIs.taskResponsiveP.tsv');

        % Require stats + region + roiTypes
        hasStat = isfile(statFileA) || isfile(statFileB) || isfile(statTsvA) || isfile(statTsvB);
        needResp = opt.responsiveOnly;

        if ~(hasStat && isfile(regFile) && isfile(typeFile) && (~needResp || isfile(respFile)))
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
            % Otherwise load NPY
            if isfile(statFileA), npyPath = statFileA; else, npyPath = statFileB; end
            S = double(readNPY(npyPath));
            if ndims(S) ~= 2
                error('Expected a 2D array in %s (nROIs x nTests).', npyPath);
            end
            testNames = "Test" + string(1:size(S,2));
        end

        ix_sel = startsWith(testNames, opt.stat);
        testNames = testNames(ix_sel);
        S = S(:, ix_sel);

        n = size(S,1);

        % -------- Load responsive p-values --------
        respNames = string.empty(1,0);
        Resp = [];
        if isfile(respFile)
            respTab = readtable(respFile, 'FileType','text', 'Delimiter','\t', 'VariableNamingRule','preserve');
            respNames = string(respTab.Properties.VariableNames);
            Resp = table2array(respTab);
            if size(Resp,1) ~= n
                warning('Size mismatch (responsiveP vs stats) in %s; ignoring responsive filter for this FOV.', fovPath);
                respNames = string.empty(1,0);
                Resp = [];
            end
        elseif opt.responsiveOnly
            warning('Missing taskResponsiveP.tsv in %s; skipping FOV because responsiveOnly=true.', fovPath);
            continue
        end

        % Update union of all tests; if new tests appear, add NaN columns to existing records
        newTests = setdiff(testNames, allTestNames, 'stable');
        if ~isempty(newTests)
            for nt = 1:numel(newTests)
                records.(newTests(nt)) = nan(height(records), 1);
                records.(newTests(nt) + "_respP") = nan(height(records), 1);
            end
            allTestNames = [allTestNames, newTests];
        end

        % Load regions + roiTypes
        regionIds = double(readNPY(regFile));
        regionIds = regionIds(:);

        roiTypes = double(readNPY(typeFile));
        roiTypes = roiTypes(:);

        if numel(regionIds) ~= n || numel(roiTypes) ~= n
            warning('Size mismatch (region/roiTypes vs stats) in %s; skipping FOV.', fovPath);
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
        T.roiType    = roiTypes;

        % Add the test columns that exist in this file
        for k = 1:numel(testNames)
            tn = testNames(k);
            T.(tn) = S(:,k);

            % matched responsive p-values for this tuning test
            respCol = tn + "_respP";
            ridx = map_tuning_test_to_responsive_test(tn, respNames);
            if ~isempty(ridx)
                T.(respCol) = Resp(:,ridx);
            else
                T.(respCol) = nan(n,1);
            end
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
                if isnumeric(T.(v))
                    records.(v) = nan(height(records), 1);
                elseif isdatetime(T.(v))
                    records.(v) = NaT(height(records), 1);
                else
                    records.(v) = repmat("", height(records), 1);
                end
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
    error('No usable FOVs found. (Need taskTunedStat.npy/tsv + brainLocation .npy + mpciROITypes.npy at minimum)');
end

% -----------------------------
% Day ordering (YYYY-MM-DD)
% -----------------------------
records.day = datetime(records.day, 'InputFormat','yyyy-MM-dd');
days = unique(records.day);
days = sort(days);
nDays = numel(days);

% Tests
metaCols = ["day","sessionID","fov","clusterUID","regionId","roiType"]; %#ok<NASGU>
testCols = allTestNames;
nTests = numel(testCols);

%% Tracked-only or all-type-1 filtering
if opt.trackedOnly
    % Original behavior: tracked ROIs only, then require chronic tracking
    hasUID = records.clusterUID ~= "" & ~ismissing(records.clusterUID);
    Rbase = records(hasUID, :);

    if isempty(Rbase)
        error('No tracked ROIs found (no non-empty clusterUIDs.csv entries in provided sessions).');
    end

    % Sanity: duplicates of same UID within the same day
    [Gdup, ~, ~] = findgroups(Rbase.day, Rbase.clusterUID);
    cnt = splitapply(@numel, Rbase.clusterUID, Gdup);
    if any(cnt > 1)
        warning('Found UID×day groups with >1 occurrence. Your chronic definition will exclude these UIDs.');
    end

    % Get chronic UIDs tracked on all days
    cUIDs_d1 = get_clusterUIDs(sessionPaths(1));
    is_tracked = get_trackedROIs(cUIDs_d1, sessionPaths(2:end));
    chronicUIDs = cUIDs_d1(is_tracked);

    Rfilt = Rbase(ismember(Rbase.clusterUID, chronicUIDs), :);

    if isempty(chronicUIDs)
        warning('No chronic UIDs found under the "exactly once per day" criterion.');
    end

    % Stable region per chronic UID = modal region across days
    chUIDs = chronicUIDs(:);
    nCh = numel(chUIDs);
    stableRegion = nan(nCh,1);
    tieFlag = false(nCh,1);

    for j = 1:nCh
        uid = chUIDs(j);
        r = Rfilt.regionId(Rfilt.clusterUID == uid);
        [stableRegion(j), tieFlag(j)] = modal_region(r);
    end

    uid2reg = table(chUIDs, stableRegion, tieFlag, ...
        'VariableNames', {'clusterUID','stableRegionId','tieFlag'});

    if nCh > 0
        [tf, loc] = ismember(Rfilt.clusterUID, uid2reg.clusterUID);
        if ~all(tf), error('UID lookup failed (unexpected).'); end
        Rfilt.stableRegionId = uid2reg.stableRegionId(loc);
    else
        Rfilt.stableRegionId = [];
    end

    recordsTracked = Rbase;
else
    % all somatic ROIs regardless of UID
    Rfilt = records(records.roiType == 1, :);
    if isempty(Rfilt)
        error('No ROIs with mpciROITypes==1 found.');
    end
    
    % Use instantaneous regionId as grouping region
    Rfilt.stableRegionId = Rfilt.regionId;

    chronicUIDs = string.empty(0,1);
    uid2reg = table(string.empty(0,1), [], false(0,1), ...
        'VariableNames', {'clusterUID','stableRegionId','tieFlag'});

    recordsTracked = table();
end

% Optional region filter
[regionIdsExpanded, regsLabel] = expand_region_ids(opt);
if ~isempty(regionIdsExpanded)
    keep = ismember(double(Rfilt.stableRegionId), double(regionIdsExpanded));
    Rfilt = Rfilt(keep, :);

    if opt.trackedOnly
        chronicUIDs = unique(Rfilt.clusterUID);
        uid2reg = uid2reg(ismember(uid2reg.clusterUID, chronicUIDs), :);
    end

    if isempty(Rfilt)
        warning('Region filter matched no ROIs. Requested: %s', regsLabel);
    end
else
    regsLabel = "All";
end

% enforce minimum N per DAY x region
minROIsPerRegion = opt.minROIsPerRegion;

% First get all candidate regions after optional region filtering
if ~isempty(Rfilt)
    regionListAll = sort(unique(double(Rfilt.stableRegionId)));
else
    regionListAll = [];
end

nRegionsAll = numel(regionListAll);

% Count ROIs per day x region
nPerDayRegion = zeros(nDays, nRegionsAll);

for di = 1:nDays
    Rd = Rfilt(Rfilt.day == days(di), :);

    for ri = 1:nRegionsAll
        rID = regionListAll(ri);
        regMask = (double(Rd.stableRegionId) == rID);

        if opt.trackedOnly
            % in tracked mode count unique tracked ROIs present that day
            u = unique(Rd.clusterUID(regMask));
            u = u(u ~= "" & ~ismissing(u));
            nPerDayRegion(di,ri) = numel(u);
        else
            % in all-ROI mode count all accepted ROIs present that day
            nPerDayRegion(di,ri) = sum(regMask);
        end
    end
end

% Keep a region if ANY day passes threshold
keepRegion = any(nPerDayRegion >= minROIsPerRegion, 1);
regionList = regionListAll(keepRegion);
nPerDayRegion = nPerDayRegion(:, keepRegion);

if isempty(regionList)
    warning('No regions have >= %d ROIs on any day after filtering.', minROIsPerRegion);
end

% Optional summary count for labels: maximum daily count in that region
if ~isempty(nPerDayRegion)
    nPerReg = max(nPerDayRegion, [], 1);
else
    nPerReg = [];
end

nRegions = numel(regionList);

%% 3) FULL DISTRIBUTIONS per day × region × test
distStats = cell(nDays, nRegions, nTests);
nDen = zeros(nDays, nRegions, nTests);

for di = 1:nDays
    Rd = Rfilt(Rfilt.day == days(di), :);

    for ri = 1:nRegions
        rID = regionList(ri);
        regMask = (Rd.stableRegionId == rID);

        % Skip this day/region unless it passes the minimum count
        if nPerDayRegion(di,ri) < minROIsPerRegion
            continue
        end

        if ~any(regMask), continue; end

        for ti = 1:nTests
            tn = testCols(ti);

            x = Rd.(tn);
            x = x(regMask);

            ok = ~isnan(x);

            % Optional: keep only task-responsive ROIs for the mapped responsive test
            if opt.responsiveOnly
                respCol = tn + "_respP";
                if ismember(respCol, string(Rd.Properties.VariableNames))
                    rp = Rd.(respCol);
                    rp = rp(regMask);
                    ok = ok & ~isnan(rp) & (rp < opt.responsiveAlpha);
                else
                    ok = false(size(ok));
                end
            end

            nDen(di,ri,ti) = sum(ok);

            if any(ok)
                distStats{di,ri,ti} = x(ok);
            else
                distStats{di,ri,ti} = [];
            end
        end
    end
end

%% 4) Plot compact per-region distribution summaries

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
for ii = 1:numel(regionList)
    idx = find(rid_all == regionList(ii), 1);
    if ~isempty(idx)
        regionLabels(ii) = acr_all(idx-1); % kept from your original
    else
        regionLabels(ii) = string(regionList(ii));
    end
end

% sort regions by structure tree row index
if nRegions > 0 && ~isempty(regionList)
    [~, idxRow] = ismember(double(regionList), rid_all);
    [~, sortIdx] = sort(idxRow);

    regionList   = regionList(sortIdx);
    regionLabels = regionLabels(sortIdx);

    distStats = distStats(:, sortIdx, :);
    nDen      = nDen(:, sortIdx, :);

    if ~isempty(nPerReg)
        nPerReg = nPerReg(sortIdx);
    end
end

for ti = 1:nTests

    % Keep only days with any data for this test
    dayHasData = false(nDays,1);
    for di = 1:nDays
        hasAny = false;
        for ri = 1:nRegions
            x = distStats{di,ri,ti};
            if ~isempty(x)
                hasAny = true;
                break
            end
        end
        dayHasData(di) = hasAny;
    end

    if isempty(opt.statRange)
        vals = distStats(:,:,ti);
        lo = min(cellfun(@(x) min([x; Inf], [], 'all'), vals), [], 'all');
        hi = max(cellfun(@(x) max([x; -Inf], [], 'all'), vals), [], 'all');
        if strcmp(opt.stat,'ccu')
            abs_hi = 1.01 * max(abs([lo hi] - 0.5));
            statRange = 0.5 + [-abs_hi, abs_hi];
        else
            abs_hi = 1.01 * max(abs([lo hi]));
            statRange = [-abs_hi, abs_hi];
        end
    else
        statRange = opt.statRange;
    end

    dayIdx = find(dayHasData);
    if isempty(dayIdx)
        warning('No data for test %s; skipping plot.', testCols(ti));
        continue
    end

    regionIdxAll = 1:nRegions;
    regionChunks = chunk_indices(regionIdxAll, opt.maxRegionsPerFigure);

    for ci = 1:numel(regionChunks)
        regionIdx = regionChunks{ci};
        nCols = numel(regionIdx);

        figure('Color','w', ...
               'Name', sprintf('%s - %s', subjStr, testCols(ti)), ...
               'Position', [60 120 max(900, 220*nCols) max(500, 28*numel(dayIdx) + 140)]);

        tl = tiledlayout(1, nCols, 'TileSpacing','compact', 'Padding','compact');

        for c = 1:nCols
            ri = regionIdx(c);
            ax = nexttile(tl, c);
            hold(ax, 'on');

            % vertical reference line
            if strcmp(opt.stat,'ccMI')
                xline(ax, 0, 'k-', 'LineWidth', 1);
            elseif strcmp(opt.stat,'ccu')
                xline(ax, 0.5, 'k-', 'LineWidth', 1);
            end

            for jj = 1:numel(dayIdx)
                di = dayIdx(jj);

                x = distStats{di,ri,ti};
                if isempty(x)
                    continue
                end

                x = x(:);
                x = x(~isnan(x) & isfinite(x));
                x = x(x >= statRange(1) & x <= statRange(2));
                if isempty(x)
                    continue
                end

                y0 = jj;

                q = quantile(x, [0.05, 0.25 0.5 0.75, 0.95]);
                outliers = x(x < q(1) | x > q(end));

                patch(ax, [q(1) q(5) q(5) q(1)], [y0-.4 y0-.4 y0+.4 y0+.4], [0.75 0.75 0.75], 'EdgeColor','none');
                patch(ax, [q(2) q(4) q(4) q(2)], [y0-.3 y0-.3 y0+.3 y0+.3], [0.35 0.35 0.35], 'EdgeColor','none');
                plot(ax, q(3), y0, 'ko', 'MarkerFaceColor', 'k');
                plot(ax, outliers, y0 * ones(size(outliers)), 'k.');
            end

            ax.YDir = 'reverse';
            ylim(ax, [0.5, numel(dayIdx) + 0.5]);
            xlim(ax, statRange);
            ax.Box = 'on';

            % Title = region
            if ~isempty(nPerReg)
                title(ax, sprintf('%s (max n=%d)', regionLabels(ri), nPerReg(ri)), ...
                    'Interpreter','none');
            else
                title(ax, regionLabels(ri), 'Interpreter','none');
            end

            % Y labels only on first subplot
            ax.YTick = 1:numel(dayIdx);
            if c == 1
                ax.YTickLabel = string(days(dayIdx), 'yyyy-MM-dd');
                ylabel(ax, 'Day');
            else
                ax.YTickLabel = {};
            end
        end

        xlabel(tl, opt.stat);
        if opt.responsiveOnly
            modeLabel = sprintf('responsiveOnly (p<%.3g)', opt.responsiveAlpha);
        else
            modeLabel = 'all somatic ROIs';
        end
        title(tl, sprintf('%s\n%s | %s', subjStr, testCols(ti), modeLabel), 'Interpreter','none');
    end
end

% Output
out = struct();
out.recordsAll      = records;
out.recordsTracked  = recordsTracked;
out.recordsFiltered = Rfilt;
out.days            = days;
out.regionList      = regionList;
out.regionLabels    = regionLabels;
out.testCols        = testCols;
out.distStats       = distStats;   % cell(nDays,nRegions,nTests) of vectors
out.nDen            = nDen;
out.chronicUIDs     = chronicUIDs;
out.uid2reg         = uid2reg;
out.nPerReg         = nPerReg;
out.nPerDayRegion   = nPerDayRegion;
out.minROIsPerRegion= minROIsPerRegion;
out.trackedOnly     = opt.trackedOnly;
out.responsiveOnly = opt.responsiveOnly;
out.responsiveAlpha = opt.responsiveAlpha;
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

%% helper: chunk indices into groups of at most maxCols
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

function draw_horizontal_violin(ax, x, y0, halfWidth, xlimRange)
% Draw a horizontal violin centered at y = y0
% x = data values
% halfWidth = max half-height of violin
% xlimRange = [xmin xmax]

x = x(:);
x = x(~isnan(x) & isfinite(x));

if isempty(x)
    return
end

xmin = xlimRange(1);
xmax = xlimRange(2);

% Strictly clip to support BEFORE ksdensity
x = x(x >= xmin & x <= xmax);

if isempty(x)
    return
end

if numel(x) == 1
    plot(ax, x, y0, 'k.', 'MarkerSize', 8);
    return
end

% If all values identical, KDE can behave badly; draw a thin marker instead
if range(x) == 0
    xv = x(1);
    plot(ax, [xv xv], [y0-halfWidth*0.6, y0+halfWidth*0.6], 'k-', 'LineWidth', 1.5);
    return
end

% evaluation grid
xi = linspace(xmin, xmax, 200);

% KDE on bounded support
[f, xi] = ksdensity(x, xi, 'Support', [xmin xmax], 'BoundaryCorrection', 'reflection');

if isempty(f) || all(~isfinite(f)) || max(f) <= 0
    plot(ax, median(x), y0, 'k.', 'MarkerSize', 8);
    return
end

f = f ./ max(f);
w = halfWidth * f;

X = [xi, fliplr(xi)];
Y = [y0 + w, fliplr(y0 - w)];

patch(ax, X, Y, [0.3 0.3 0.3], ...
    'FaceAlpha', 0.35, ...
    'EdgeColor', 'none');
end

