function out = mpci_tuning_heatmaps(sessionPaths, alpha, opt)
% mpci_tuning_heatmaps plots task-variable tuning heatmaps divided by brain
% region and session
%
%   %TODO update header 
%   - Loads p-values from mpciROIs.taskTunedP.tsv
%   - Computes fraction significant (two-tailed) per day × region × test
%   - Plots heatmap (dates on X, regions on Y) per test
%   - Also loads test statistics from mpciROIs.taskTunedStat.tsv
%   - Computes MEAN statistic per day × region × test (across chronic ROIs)
%   - Adds a subplot underneath the fraction-significant heatmap showing
%     the mean statistic heatmap with a diverging RdBu-like colormap centered at 0.

if nargin < 2 || isempty(alpha), alpha = 0.05; end
if nargin < 3, opt = struct(); end

% optional region filter inputs
if ~isfield(opt,'region'), opt.region =  {'VIS*','RSP*','SS*','AUD*','MO*'}; end % acronyms/tokens
if ~isfield(opt,'st'),     opt.st = ''; end
if ~isfield(opt,'stPath'), opt.stPath = "C:\Users\Samuel\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv"; end
if ~isfield(opt,'trackedOnly') || isempty(opt.trackedOnly), opt.trackedOnly = true; end

% optional: fixed symmetric color range for the STAT mean heatmap (else auto per test)
% e.g. opt.statCMax = 0.5;  -> caxis([-0.5 0.5])
if ~isfield(opt,'statCMax'), opt.statCMax = []; end

sessionPaths = string(sessionPaths(:));
if isempty(sessionPaths)
    error('sessionPaths is empty.');
end

records = table();
allTestNames = string.empty(1,0);   % union of all test columns encountered (from P TSV)

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

        pFile    = fullfile(fovPath, 'mpciROIs.taskTunedP.tsv');
        statFile = fullfile(fovPath, 'mpciROIs.taskTunedStat.tsv');
        regFile  = fullfile(fovPath, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy');
        uidFile  = fullfile(fovPath, 'mpciROIs.clusterUIDs.csv');
        typeFile = fullfile(fovPath, 'mpciROIs.mpciROITypes.npy');

        % Always require p + region. Stats is optional (for the added subplot).
        if ~(isfile(pFile) && isfile(regFile) && isfile(typeFile))
            continue
        end

        % Load P-values
        pTab = readtable(pFile, 'FileType','text', 'Delimiter','\t', 'VariableNamingRule','preserve');
        testNames = string(pTab.Properties.VariableNames);
        P = table2array(pTab);
        n = size(P,1);

        % Load STAT table if present
        haveStat = isfile(statFile);
        if haveStat
            sTab = readtable(statFile, 'FileType','text', 'Delimiter','\t', 'VariableNamingRule','preserve');
            statNames = string(sTab.Properties.VariableNames);
            S = table2array(sTab);
            if size(S,1) ~= n
                warning('Size mismatch (stats vs P) in %s; ignoring stats for this FOV.', fovPath);
                haveStat = false;
                statNames = string.empty(1,0);
                S = [];
            end
        else
            statNames = string.empty(1,0);
            S = [];
        end

        % Update union of all tests; if new tests appear, add NaN columns to existing records
        newTests = setdiff(testNames, allTestNames, 'stable');
        if ~isempty(newTests)
            for nt = 1:numel(newTests)
                records.(newTests(nt)) = nan(height(records), 1);
                records.(newTests(nt) + "_stat") = nan(height(records), 1); % NEW: paired stat column
            end
            allTestNames = [allTestNames, newTests];
        end

        % Load roiTypes
        roiTypes = double(readNPY(typeFile));
        roiTypes = roiTypes(:);

        % Load regions
        regionIds = double(readNPY(regFile));
        regionIds = regionIds(:);

        if numel(regionIds) ~= n || numel(roiTypes) ~= n
            warning('Size mismatch (region/roiTypes vs P) in %s; skipping FOV.', fovPath);
            continue
        end

        % Load UIDs if present; otherwise fill empty strings
        if isfile(uidFile)
            uidCol = read_uid_csv(uidFile);
            if numel(uidCol) ~= n
                warning('Size mismatch (UID vs P) in %s; filling empty UIDs for this FOV.', fovPath);
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

        % Add test columns that exist in this P file
        % ALSO: add paired "_stat" columns (aligned by matching variable names when possible)
        for k = 1:numel(testNames)
            tn = testNames(k);
            T.(tn) = P(:,k);

            % stats: if present and has matching column name, use it; else NaN
            statColName = tn + "_stat";
            if haveStat
                idxS = find(statNames == tn, 1);
                if ~isempty(idxS)
                    T.(statColName) = S(:,idxS);
                else
                    T.(statColName) = nan(n,1);
                end
            else
                T.(statColName) = nan(n,1);
            end
        end

        % ---- Robust append to records with variable columns ----
        if isempty(records)
            records = T;
            allTestNames = testNames;
        else
            allTestNames = unique([allTestNames, testNames], 'stable');

            % 1) Add any columns present in T but missing in records
            missingInRecords = setdiff(T.Properties.VariableNames, records.Properties.VariableNames, 'stable');
            for m = 1:numel(missingInRecords)
                v = missingInRecords{m};
                if endsWith(v, "_stat")
                    records.(v) = nan(height(records), 1);
                else
                    % could be meta columns or test columns
                    if ismember(string(v), ["day","sessionID","fov","clusterUID","regionId"])
                        records.(v) = repmat("", height(records), 1); % will be overwritten by type fix below if needed
                    else
                        records.(v) = nan(height(records), 1);
                    end
                end
            end

            % 2) Add any columns present in records but missing in T
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
    error('No usable FOVs found. (Need taskTunedP.tsv + brainLocation .npy at minimum)');
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

if opt.trackedOnly
    % tracked ROIs only, then require chronic tracking
    hasUID = records.clusterUID ~= "" & ~ismissing(records.clusterUID);
    Rbase = records(hasUID, :);

    if isempty(Rbase)
        error('No tracked ROIs found (no non-empty clusterUIDs.csv entries in provided sessions).');
    end

    % Sanity: duplicates of same UID within the same day
    [G, ~, ~] = findgroups(Rbase.day, Rbase.clusterUID);
    cnt = splitapply(@numel, Rbase.clusterUID, G);
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

else
    % all somatic ROIs regardless of UID
    Rfilt = records(records.roiType == 1, :);

    if isempty(Rfilt)
        error('No ROIs with mpciROITypes==1 found.');
    end

    % Use instantaneous regionId as the grouping region
    Rfilt.stableRegionId = Rfilt.regionId;

    chronicUIDs = string.empty(0,1);
    uid2reg = table(string.empty(0,1), [], false(0,1), ...
        'VariableNames', {'clusterUID','stableRegionId','tieFlag'});
end

%% 2) Stable region per chronic UID = modal region across days
chUIDs = chronicUIDs(:);
nCh = numel(chUIDs);
stableRegion = nan(nCh,1);
tieFlag = false(nCh,1);

for j = 1:nCh
    uid = chUIDs(j);
    r = Rfilt.regionId(Rfilt.clusterUID == uid); % should be nDays long
    [stableRegion(j), tieFlag(j)] = modal_region(r);
end

uid2reg = table(chUIDs, stableRegion, tieFlag, ...
    'VariableNames', {'clusterUID','stableRegionId','tieFlag'});

if nCh > 0
    [tf, loc] = ismember(Rfilt.clusterUID, uid2reg.clusterUID);
    if ~all(tf), error('UID lookup failed (unexpected).'); end
    Rfilt.stableRegionId = uid2reg.stableRegionId(loc);

    % Optional region filter (by stable modal region per UID)
    [regionIdsExpanded, regsLabel] = expand_region_ids(opt);

    if ~isempty(regionIdsExpanded)
        keep = ismember(double(Rfilt.stableRegionId), double(regionIdsExpanded));
        Rfilt = Rfilt(keep, :);

        chronicUIDs = unique(Rfilt.clusterUID);
        uid2reg = uid2reg(ismember(uid2reg.clusterUID, chronicUIDs), :);

        if isempty(Rfilt)
            warning('Region filter matched no chronic UIDs. Requested: %s', regsLabel);
        end
    else
        regsLabel = "All";
    end
else
    Rfilt.stableRegionId = [];
    regsLabel = "All";
end

% enforce minimum N ROIs per region
minROIsPerRegion = 50;
nUIDsPerReg = [];

if ~isempty(Rfilt)
    [G, regVals] = findgroups(double(Rfilt.stableRegionId));

    if opt.trackedOnly
        nPerReg = splitapply(@(u) numel(unique(u)), Rfilt.clusterUID, G);
    else
        nPerReg = splitapply(@numel, Rfilt.regionId, G);
    end

    keepReg = nPerReg >= minROIsPerRegion;
    keepRegionIds = regVals(keepReg);

    Rfilt = Rfilt(ismember(double(Rfilt.stableRegionId), keepRegionIds), :);

    nPerReg = nPerReg(keepReg);
    keepRegionIds = regVals(keepReg);

    % Filter Rfilt to keep only those regions
    Rfilt = Rfilt(ismember(double(Rfilt.stableRegionId), keepRegionIds), :);

    % Keep these consistent downstream (optional but recommended)
    chronicUIDs = unique(Rfilt.clusterUID);
    uid2reg = uid2reg(ismember(uid2reg.clusterUID, chronicUIDs), :);
    nUIDsPerReg = nUIDsPerReg(keepReg);

    if isempty(keepRegionIds)
        warning('No regions have >= %d ROIs after filtering.', minROIsPerRegion);
    end
end

% Regions based on stable IDs (after filtering)
if ~isempty(Rfilt)
    regionList = sort(unique(Rfilt.stableRegionId));
else
    regionList = [];
end
nRegions = numel(regionList);

%% 3) Fractions significant per day × stable region × test
% AND mean statistic per day × stable region × test (from *_stat columns)

frac = nan(nDays, nRegions, nTests);
meanStat = nan(nDays, nRegions, nTests);

nDen = zeros(nDays, nRegions, nTests);       % denom for P-values
nDenStat = zeros(nDays, nRegions, nTests);   % denom for stats (may differ if NaNs)

for di = 1:nDays
    if nCh == 0, continue; end
    Rd = Rfilt(Rfilt.day == days(di), :);

    for ri = 1:nRegions
        rID = regionList(ri);
        regMask = (Rd.stableRegionId == rID);

        if ~any(regMask), continue; end

        for ti = 1:nTests
            tn = testCols(ti);

            % P-values
            p = Rd.(tn);
            p = p(regMask);

            okp = ~isnan(p);
            denomP = sum(okp);
            nDen(di,ri,ti) = denomP;

            if denomP > 0
                sig = (p(okp) < alpha/2) | (p(okp) > 1 - alpha/2);
                frac(di,ri,ti) = sum(sig) / denomP;
            end

            % Stats mean (paired column tn+"_stat")
            sn = tn + "_stat";
            if ismember(sn, string(Rd.Properties.VariableNames))
                s = Rd.(sn);
                s = s(regMask);
                oks = ~isnan(s);
                nDenStat(di,ri,ti) = sum(oks);
                if any(oks)
                    meanStat(di,ri,ti) = mean(s(oks));
                end
            end
        end
    end
end

%% 4) Plot heatmaps (top: fraction significant, bottom: mean stat centered at 0)

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

for i = 1:numel(regionList)
    idx = find(rid_all == regionList(i), 1);
    region_ix(i) = idx;
    if ~isempty(idx)
        regionLabels(i) = acr_all(idx-1); % hack to remove '1' (kept as-is)
    else
        regionLabels(i) = string(regionList(i));
    end
end

% sorting of region (by allen structure tree row index)
if nRegions > 0
    [tf, idxRow] = ismember(double(regionList), rid_all); %#ok<ASGLU>
    [~, sortIdx] = sort(idxRow);
    regionLabels = regionLabels(sortIdx);
    regionList = regionList(sortIdx);
    frac = frac(:, sortIdx, :);
    meanStat = meanStat(:, sortIdx, :);
    nDen = nDen(:, sortIdx, :);
    nDenStat = nDenStat(:, sortIdx, :);
    if ~isempty(nPerReg)
        nPerReg = nPerReg(sortIdx);
    end
end

for ti = 1:nTests
    % Keep same "dayIdx" logic as original (based on P-values)
    dayHasData = squeeze(sum(nDen(:,:,ti), 2)) > 0;  % nDays x 1 (sum across regions)
    dayIdx = find(dayHasData);

    if isempty(dayIdx)
        warning('No data for test %s; skipping plot.', testCols(ti));
        continue
    end

    Mfrac = frac(dayIdx,:,ti).';      % (nRegions x nDays_for_this_test)
    Mstat = meanStat(dayIdx,:,ti).';  % (nRegions x nDays_for_this_test)

    % symmetric caxis centered at 0
    if ~isempty(opt.statCMax)
        cmax = abs(opt.statCMax);
    else
        cmax = max(abs(Mstat(:)), [], 'omitnan');
        if isempty(cmax) || ~isfinite(cmax) || cmax == 0
            cmax = 1; % fallback
        end
    end

    figure('Color','w', 'Name', char(testCols(ti)), ...
        'Position',[681,200,560, max(260, nRegions*24 + 200)]);

    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    % ---- TOP: fraction significant heatmap ----
    ax1 = nexttile(tl, 1);
    imagesc(ax1, Mfrac);
    colormap(ax1, parula);
    cb1 = colorbar(ax1); %#ok<NASGU>

    ax1.XTick = 1:numel(dayIdx);
    ax1.XTickLabel = string(days(dayIdx), 'yyyy-MM-dd');
    ax1.XTickLabelRotation = 45;

    i_firstBias = find(days(dayIdx) >= getSubjectDate(subjStr), 1, 'first');
    if ~isempty(i_firstBias)
        xline(ax1, i_firstBias-0.5,'r--','linewidth',2);
    end

    ax1.YTick = 1:max(1,nRegions);
    if nRegions > 0
        ax1.YTick = 1:nRegions;
        ax1.YTickLabel = regionLabels + " (" + string(nUIDsPerReg) + ")";
    else
        ax1.YTickLabel = {};
    end

    ylabel(ax1, 'Region');
    title(ax1, sprintf('%s - Fraction sign. (alpha=%.3g)\n%s', ...
        subjStr, alpha, testCols(ti)), 'Interpreter','none');

    % ---- BOTTOM: mean statistic heatmap ----
    ax2 = nexttile(tl, 2);
    imagesc(ax2, Mstat);
    if strncmpi(testCols(ti),'ccMI',4)
        colormap(ax2, brewermap(256,'*RdBu'));
        caxis(ax2, [-cmax, cmax]); 
    elseif strncmpi(testCols(ti),'ccu',3)
        colormap(ax2, brewermap(256,'*RdBu'));
        cmax = 0.05;
        %cmax = max(abs(Mstat(:)-0.5), [], 'omitnan');
        caxis(ax2, 0.5+[-cmax, cmax])
    end
    cb2 = colorbar(ax2); %#ok<NASGU>

    ax2.XTick = 1:numel(dayIdx);
    ax2.XTickLabel = string(days(dayIdx), 'yyyy-MM-dd');
    ax2.XTickLabelRotation = 45;

    if ~isempty(i_firstBias)
        xline(ax2, i_firstBias-0.5,'r--','linewidth',2);
    end

    ax2.YTick = 1:max(1,nRegions);
    if nRegions > 0
        ax2.YTick = 1:nRegions;
        ax2.YTickLabel = regionLabels + " (" + string(nPerReg) + ")";
    else
        ax2.YTickLabel = {};
    end

    xlabel(ax2, 'Day');
    ylabel(ax2, 'Region');
    title(ax2, sprintf('%s - Mean test statistic\n%s', ...
        subjStr, testCols(ti)), 'Interpreter','none');

end

% Output
out = struct();
out.recordsAll = records;
out.recordsFiltered = Rfilt;
out.trackedOnly = opt.trackedOnly;
out.nPerReg = nPerReg;
out.days = days;
out.regionList = regionList;
out.testCols = testCols;
out.frac = frac;
out.meanStat = meanStat;         % NEW
out.nDen = nDen;
out.nDenStat = nDenStat;         % NEW
out.chronicUIDs = chronicUIDs;
out.uid2reg = uid2reg;
out.alpha = alpha;

end

%% helper
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
