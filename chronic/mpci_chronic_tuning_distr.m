function out = mpci_chronic_tuning_distr(sessionPaths, opt)
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
if ~isfield(opt,'statRange') || isempty(opt.statRange), opt.statRange = []; end
if ~isfield(opt,'violinWidth') || isempty(opt.violinWidth), opt.violinWidth = .38; end  % optional: cap cols shown at once
if ~isfield(opt,'stat') || isempty(opt.stat), opt.stat = 'ccMI'; end
if ~isfield(opt,'plotStyle') || isempty(opt.plotStyle), opt.plotStyle = 'box'; end
if ~isfield(opt,'cv') || isempty(opt.cv), opt.cv = false; end
% valid: 'violin', 'box', 'quantileLines'

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

        ix_sel = startsWith(testNames,opt.stat);

        testNames = testNames(ix_sel);
        S = S(:,ix_sel);

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
allCols = allTestNames(:)';
isBase = ~endsWith(allCols, "_odd") & ~endsWith(allCols, "_even");
testCols = allCols(isBase);
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

%% 4) Plot per-region distribution summaries
% plot styles:
%   'violin'        : one subplot per region, one violin per day
%   'box'           : one subplot per region, one box-summary per day
%   'quantileLines' : one tall subplot per region, days on x-axis, one line per quantile

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

    if ~isempty(nUIDsPerReg)
        nUIDsPerReg = nUIDsPerReg(sortIdx);
    end
end

% Optional display settings
if ~isfield(opt,'maxRegionsPerFigure') || isempty(opt.maxRegionsPerFigure)
    opt.maxRegionsPerFigure = 10;   % paginate if many regions
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
        if strcmp(opt.plotStyle,'quantileLines')
            lo = min(cellfun(@(x) quantile([x; Inf],0.1, 'all'), distStats(:,:,ti)), [], 'all');
            hi = max(cellfun(@(x) quantile([x; -Inf],0.9, 'all'), distStats(:,:,ti)), [], 'all');
        else
            lo = min(cellfun(@(x) min([x; Inf], [], 'all'), distStats(:,:,ti)), [], 'all');
            hi = max(cellfun(@(x) max([x; -Inf], [], 'all'), distStats(:,:,ti)), [], 'all');
        end
        if strcmp(opt.stat,'ccu')
            abs_hi = 1.01*max(abs([lo hi]-0.5));
            statRange = 0.5+[-abs_hi, abs_hi];
        else
            abs_hi = 1.01*max(abs([lo hi]));
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

    switch lower(string(opt.plotStyle))

        case "quantilelines"
            % One tall figure with vertically stacked subplots (one per region)
            %
            % If opt.cv == false:
            %   - define quantile groups from the reference day's base test values
            %   - follow the same ROIs across days
            %   - plot the median of each group across days
            %
            % If opt.cv == true:
            %   - define quantile groups from the reference day's ODD values
            %   - follow the same ROIs across days
            %   - plot the median of each group across days using EVEN values
            %   - if odd/even columns are missing, skip this test

            regionIdxAll = 1:nRegions;
            regionChunks = chunk_indices(regionIdxAll, opt.maxRegionsPerFigure);

            % define ref-day quantiles
            qProb = [0 0.1:0.1:0.9 1];
            nGroups = numel(qProb)-1;
            lineColors = brewermap(nGroups, '*RdBu');

            biasDay = getSubjectDate(subjStr);

            tnBase = testCols(ti);

            if opt.cv
                tnRef  = tnBase + "_odd";
                tnPlot = tnBase + "_even";

                allVarNames = string(Rc.Properties.VariableNames);
                if ~ismember(tnRef, allVarNames) || ~ismember(tnPlot, allVarNames)
                    warning('Skipping %s: opt.cv=true but missing %s and/or %s.', tnBase, tnRef, tnPlot);
                    continue
                end
            else
                tnRef  = tnBase;
                tnPlot = tnBase;
            end

            % Compute valid days separately for this test
            validDay = false(numel(dayIdx), 1);
            for jj = 1:numel(dayIdx)
                di = dayIdx(jj);
                Rd = Rc(Rc.day == days(di), :);
                if isempty(Rd), continue; end

                xRefDay = Rd.(tnRef);

                if opt.cv
                    xPlotDay = Rd.(tnPlot);
                    validDay(jj) = any(~isnan(xRefDay)) && any(~isnan(xPlotDay));
                else
                    validDay(jj) = any(~isnan(xRefDay));
                end
            end

            % First valid day on/after biasDay for this specific test
            if true
                validOnOrAfterBias = validDay & (days(dayIdx) >= biasDay);
                biasOrd = find(validOnOrAfterBias, 1, 'first');
                if isempty(biasOrd)
                    warning('No valid reference day found on/after bias day for test %s; skipping.', tnBase);
                    continue
                end
                diRef = dayIdx(biasOrd);
            else
                dayIdxValid = dayIdx(validDay); diRef = dayIdxValid(1);
            end


            % Plot only valid days for this test
            dayIdxPlot = dayIdx(validDay);

            if isempty(dayIdxPlot) || length(dayIdxPlot)<=1
                warning('Need > 1 valid plottin day for test %s; skipping.', tnBase);
                continue
            end

            for ci = 1:numel(regionChunks)
                regionIdx = regionChunks{ci};
                nRows = numel(regionIdx);

                figure('Color','w', ...
                    'Name', sprintf('%s - %s', subjStr, tnBase), ...
                    'Position', [120 60 300 max(350, 180*nRows + 80)]);

                tl = tiledlayout(nRows, 1, 'TileSpacing','compact', 'Padding','compact');

                for rr = 1:nRows
                    ri = regionIdx(rr);
                    ax = nexttile(tl, rr);
                    hold(ax, 'on');

                    rID = regionList(ri);

                    % Rows for this region on reference day
                    RdRef = Rc(Rc.day == days(diRef) & Rc.stableRegionId == rID, :);

                    xRef = RdRef.(tnRef);
                    uidRef = RdRef.clusterUID;

                    okRef = ~isnan(xRef) & uidRef ~= "" & ~ismissing(uidRef);
                    xRef = xRef(okRef);
                    uidRef = uidRef(okRef);

                    % clip to plotting range for consistency
                    xRef = xRef(:);
                    uidRef = uidRef(:);
                    keepRange = xRef >= statRange(1) & xRef <= statRange(2);
                    xRef = xRef(keepRange);
                    uidRef = uidRef(keepRange);

                    if numel(xRef) < nGroups
                        title(ax, sprintf('%s (too few ROIs on ref day)', regionLabels(ri)), ...
                            'Interpreter','none');
                        continue
                    end

                    % Quantile edges from reference-day distribution
                    qVals = quantile(xRef, qProb);

                    % Make bin edges robust to repeated quantiles
                    grpRef = nan(size(xRef));
                    for gi = 1:nGroups
                        lo = qVals(gi);
                        hi = qVals(gi+1);
                        if gi < nGroups
                            grpRef(xRef >= lo & xRef < hi) = gi;
                        else
                            grpRef(xRef >= lo & xRef <= hi) = gi;
                        end
                    end

                    % Fallback for repeated-edge cases: force assignment by rank
                    if any(isnan(grpRef))
                        [~, sortIx] = sort(xRef, 'ascend');
                        grpRef = nan(size(xRef));
                        nRef = numel(xRef);
                        for k = 1:nRef
                            gi = min(nGroups, ceil(k / nRef * nGroups));
                            grpRef(sortIx(k)) = gi;
                        end
                    end

                    % Store UID sets for each group
                    groupUIDs = cell(1, nGroups);
                    for gi = 1:nGroups
                        groupUIDs{gi} = uidRef(grpRef == gi);
                    end

                    % Median trajectory for each group across days
                    M = nan(nGroups, numel(dayIdxPlot));

                    for jj = 1:numel(dayIdxPlot)
                        di = dayIdxPlot(jj);
                        Rd = Rc(Rc.day == days(di) & Rc.stableRegionId == rID, :);
                        if isempty(Rd), continue; end

                        xDay = Rd.(tnPlot);
                        uidDay = Rd.clusterUID;

                        okDay = ~isnan(xDay) & uidDay ~= "" & ~ismissing(uidDay);
                        xDay = xDay(okDay);
                        uidDay = uidDay(okDay);

                        if isempty(xDay), continue; end

                        for gi = 1:nGroups
                            if isempty(groupUIDs{gi}), continue; end
                            keep = ismember(uidDay, groupUIDs{gi});
                            if any(keep)
                                M(gi,jj) = median(xDay(keep), 'omitnan');
                            end
                        end
                    end

                    % reference line
                    if strcmp(opt.stat,'ccMI')
                        yline(ax, 0, 'k-', 'LineWidth', 1);
                    elseif strcmp(opt.stat,'ccu')
                        yline(ax, 0.5, 'k-', 'LineWidth', 1);
                    end

                    % plot group trajectories
                    for gi = 1:nGroups
                        plot(ax, 1:numel(dayIdxPlot), M(gi,:), '-', ...
                            'Color', lineColors(gi,:), 'LineWidth', 1.8);
                    end

                    xlim(ax, [1 numel(dayIdxPlot)]);
                    ylim(ax, statRange);
                    ax.Box = 'on';

                    if ~isempty(nUIDsPerReg)
                        title(ax, sprintf('%s (%d)', regionLabels(ri), nUIDsPerReg(ri)), ...
                            'Interpreter','none');
                    else
                        title(ax, regionLabels(ri), 'Interpreter','none');
                    end

                    ax.XTick = 1:numel(dayIdxPlot);
                    if rr == nRows
                        ax.XTickLabel = string(days(dayIdxPlot), 'yyyy-MM-dd');
                        ax.XTickLabelRotation = 45;
                        xlabel(ax, 'Day');
                    else
                        ax.XTickLabel = {};
                    end
                    ylabel(ax, opt.stat);
                end

                if opt.cv
                    title(tl, sprintf('%s\n%s\n(groups: odd, plotted: even)', subjStr, tnBase), ...
                        'Interpreter','none');
                else
                    title(tl, sprintf('%s\n%s', subjStr, tnBase), ...
                        'Interpreter','none');
                end
            end

        otherwise
            % violin or box: original layout, one row of region subplots
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

                        switch lower(string(opt.plotStyle))
                            case "violin"
                                draw_horizontal_violin(ax, x, y0, opt.violinWidth, statRange);
                                q = quantile(x, [0.25 0.5 0.75]);
                                plot(ax, [q(1) q(3)], [y0 y0], 'k-', 'LineWidth', 2);
                                plot(ax, q(2), y0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

                            otherwise % 'box'
                                q = quantile(x, [0.05, 0.25 0.5 0.75, 0.95]);
                                outliers = x(x < q(1) | x > q(end));

                                patch(ax, [q(1) q(5) q(5) q(1)], [y0-.4 y0-.4 y0+.4 y0+.4], ...
                                    [0.75 0.75 0.75], 'EdgeColor','none');
                                patch(ax, [q(2) q(4) q(4) q(2)], [y0-.3 y0-.3 y0+.3 y0+.3], ...
                                    [0.35 0.35 0.35], 'EdgeColor','none');
                                plot(ax, q(3), y0, 'ko', 'MarkerFaceColor', 'k');
                                plot(ax, outliers, y0 * ones(size(outliers)), 'k.');
                        end
                    end

                    ax.YDir = 'reverse';
                    ylim(ax, [0.5, numel(dayIdx) + 0.5]);
                    xlim(ax, statRange);
                    ax.Box = 'on';

                    % Title = region
                    if ~isempty(nUIDsPerReg)
                        title(ax, sprintf('%s (%d)', regionLabels(ri), nUIDsPerReg(ri)), ...
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
                title(tl, sprintf('%s\n%s', subjStr, testCols(ti)), 'Interpreter','none');
            end
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
%out.histEdges       = edges;
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