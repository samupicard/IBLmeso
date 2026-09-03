function Tsess = mpci_tuning_StatsByRegion(sessionPath, alpha, opt)
% mpci_tuning_StatByRegion
% Compute per-session per-region per-test summaries using ALL ROIs
% filtered by mpciROITypes==1 (no chronic filtering).
%
% Output Tsess is LONG format with rows:
%   subjStr, day, dRelBias, sessionID, fov(="ALL"), regionId, testName,
%   fracSig, nValid, nSig, nROIs,
%   medianStat, iqrStat, madStat, nValidStat
%
% Notes:
% - Variable test set across FOVs is supported; missing tests => NaN.
% - Two-tailed criterion: p < alpha/2 OR p > 1-alpha/2.
% - Region filtering via expand_region_ids(opt) applies to regionId.
% - taskTunedStat.tsv is assumed to have the same column names/shapes as
%   taskTunedP.tsv, but matching is done by column name for robustness.
%
% Requirements: readNPY, loadStructureTree, expand_region_ids, getSubjectDate (or opt.biasDateFcn)

if nargin < 2 || isempty(alpha), alpha = 0.05; end
if nargin < 3, opt = struct(); end

% Optional options
if ~isfield(opt,'region'), opt.region = {}; end
if ~isfield(opt,'st'),     opt.st = ''; end
if ~isfield(opt,'stPath'), opt.stPath = ""; end
if ~isfield(opt,'biasDateFcn') || isempty(opt.biasDateFcn)
    opt.biasDateFcn = @getSubjectDate; % expects getSubjectDate(subjStr)
end
if ~isfield(opt,'responsiveOnly') || isempty(opt.responsiveOnly), opt.responsiveOnly = false; end
if ~isfield(opt,'responsiveAlpha') || isempty(opt.responsiveAlpha), opt.responsiveAlpha = 0.005; end

sessionPath = string(sessionPath);

% Resolve alf path
if endsWith(sessionPath, filesep + "alf") || endsWith(sessionPath, "/alf") || endsWith(sessionPath, "\alf")
    alfPath = sessionPath;
    sessRoot = fileparts(alfPath);
else
    sessRoot = sessionPath;
    alfPath = fullfile(sessRoot, 'alf');
end
if ~isfolder(alfPath)
    error('alf folder not found: %s', alfPath);
end

% Parse subj/day/sessionID from .../<subj>/<YYYY-MM-DD>/<sessionID>
[parentDir, sessionID] = fileparts(sessRoot);
[subjDir, dayStr] = fileparts(parentDir);
[~, subjStr] = fileparts(subjDir);

dayStr = string(dayStr);
subjStr = string(subjStr);
sessionID = string(sessionID);

day = datetime(dayStr, 'InputFormat','yyyy-MM-dd');
biasDay = opt.biasDateFcn(subjStr);
dRelBias = days(day - biasDay);

% Enumerate FOVs
fovs = dir(fullfile(alfPath, 'FOV*'));
fovs = fovs([fovs.isdir]);

if isempty(fovs)
    warning('No FOVs found in %s', alfPath);
    Tsess = table();
    return
end

% Accumulate per-ROI rows across FOVs with union of test columns
R = table();                 % per-ROI table (regionId + p/stat columns)
allTestNames = string.empty(1,0);

for f = 1:numel(fovs)
    fovPath = fullfile(fovs(f).folder, fovs(f).name);

    pFile    = fullfile(fovPath, 'mpciROIs.taskTunedP.tsv');
    statFile = fullfile(fovPath, 'mpciROIs.taskTunedStat.tsv');
    regFile  = fullfile(fovPath, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy');
    typeFile = fullfile(fovPath, 'mpciROIs.mpciROITypes.npy');
    respFile = fullfile(fovPath, 'mpciROIs.taskResponsiveP.tsv');

    if ~(isfile(pFile) && isfile(regFile) && isfile(typeFile) && (~opt.responsiveOnly || isfile(respFile)))
        continue
    end

    % Load P-values
    pTab = readtable(pFile, 'FileType','text', 'Delimiter','\t', 'VariableNamingRule','preserve');
    testNamesP = string(pTab.Properties.VariableNames);
    P = table2array(pTab);
    n = size(P,1);

    % Load stats table if present
    haveStat = isfile(statFile);
    if haveStat
        statTab = readtable(statFile, 'FileType','text', 'Delimiter','\t', 'VariableNamingRule','preserve');
        testNamesStat = string(statTab.Properties.VariableNames);
        S = table2array(statTab);

        if size(S,1) ~= n
            warning('Size mismatch between P and Stat tables in %s; ignoring stats for this FOV.', fovPath);
            haveStat = false;
            testNamesStat = string.empty(1,0);
            S = [];
        end
    else
        testNamesStat = string.empty(1,0);
        S = [];
    end

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

    % Union tests based on P table (primary source)
    newTests = setdiff(testNamesP, allTestNames, 'stable');
    if ~isempty(newTests)
        % add NaN columns to already accumulated R
        for nt = 1:numel(newTests)
            tn = newTests(nt);
            R.(tn) = nan(height(R), 1);              % p-value
            R.(tn + "_stat") = nan(height(R), 1);    % statistic
            R.(tn + "_respP") = nan(height(R), 1);   % responsive P
        end
        allTestNames = [allTestNames, newTests];
    end

    % Load region + roiTypes
    regionIds = double(readNPY(regFile)); regionIds = regionIds(:);
    roiTypes  = double(readNPY(typeFile)); roiTypes = roiTypes(:);

    if numel(regionIds) ~= n || numel(roiTypes) ~= n
        warning('Size mismatch in %s; skipping FOV.', fovPath);
        continue
    end

    keep = (roiTypes == 1);
    if ~any(keep), continue; end

    % Per-ROI rows for this FOV
    T = table();
    T.regionId = double(regionIds(keep));

    % Initialize all known columns
    for k = 1:numel(allTestNames)
        tn = allTestNames(k);
        T.(tn) = nan(height(T), 1);
        T.(tn + "_stat") = nan(height(T), 1);
    end

    % Fill p-value columns
    for k = 1:numel(testNamesP)
        tn = testNamesP(k);
        tmp = P(:,k);
        T.(tn) = tmp(keep);
        % matched responsive p-values for this tuning test
        respCol = tn + "_respP";
        ridx = map_tuning_test_to_responsive_test(tn, respNames);
        if ~isempty(ridx)
            T.(respCol) = Resp(keep,ridx);
        else
            T.(respCol) = nan(height(T),1);
        end
    end

    % Fill stat columns by name match
    if haveStat
        for k = 1:numel(testNamesP)
            tn = testNamesP(k);
            idxs = find(testNamesStat == tn, 1);
            if ~isempty(idxs)
                tmp = S(:, idxs);
                T.(tn + "_stat") = tmp(keep);
            end
        end
    end

    % Robust append union columns (align T and R)
    if isempty(R)
        R = T;
    else
        missingInR = setdiff(T.Properties.VariableNames, R.Properties.VariableNames, 'stable');
        for m = 1:numel(missingInR)
            v = missingInR{m};
            if isnumeric(T.(v))
                R.(v) = nan(height(R), 1);
            else
                R.(v) = repmat("", height(R), 1);
            end
        end

        missingInT = setdiff(R.Properties.VariableNames, T.Properties.VariableNames, 'stable');
        for m = 1:numel(missingInT)
            v = missingInT{m};
            if isnumeric(R.(v))
                T.(v) = nan(height(T), 1);
            else
                T.(v) = repmat("", height(T), 1);
            end
        end

        T = T(:, R.Properties.VariableNames);
        R = [R; T]; %#ok<AGROW>
    end
end

if isempty(R)
    Tsess = table();
    return
end

% Optional region filter on regionId
[regionIdsExpanded, ~] = expand_region_ids(opt);
if ~isempty(regionIdsExpanded)
    R = R(ismember(double(R.regionId), double(regionIdsExpanded)), :);
    if isempty(R)
        Tsess = table();
        return
    end
end

% Determine test columns from P columns only
allVars = string(R.Properties.VariableNames);
isStatCol = endsWith(allVars, "_stat");
isRespPCol = endsWith(allVars, "_respP");
testCols = setdiff(allVars(~isStatCol & ~isRespPCol), "regionId", 'stable');

% Compute per-region per-test summaries (long format)
regionList = unique(double(R.regionId));
nRegions = numel(regionList);
nTests = numel(testCols);

rows = nRegions * nTests;

subjCol       = repmat(subjStr, rows, 1);
dayCol        = repmat(day, rows, 1);
dRelCol       = repmat(dRelBias, rows, 1);
sessIDCol     = repmat(sessionID, rows, 1);
fovCol        = repmat("ALL", rows, 1);
regionCol     = nan(rows, 1);
testNameCol   = strings(rows, 1);

fracCol       = nan(rows, 1);
nValidCol     = zeros(rows, 1);
nSigCol       = zeros(rows, 1);
nROIsCol      = zeros(rows, 1);

medianStatCol = nan(rows, 1);
iqrStatCol    = nan(rows, 1);
madStatCol    = nan(rows, 1);
nValidStatCol = zeros(rows, 1);

rr = 0;
for ri = 1:nRegions
    rID = regionList(ri);
    regMask = (double(R.regionId) == rID);
    nROIs = sum(regMask);

    for ti = 1:nTests
        rr = rr + 1;
        tn = testCols(ti);

        % --- Fraction significant from p-values ---
        p = R.(tn);
        p = p(regMask);

        ok = ~isnan(p);

        % Optional: keep only task-responsive ROIs for the mapped responsive test
        if opt.responsiveOnly
            respCol = tn + "_respP";
            if ismember(respCol, string(R.Properties.VariableNames))
                rp = R.(respCol);
                rp = rp(regMask);
                ok = ok & ~isnan(rp) & (rp < opt.responsiveAlpha);
            else
                ok = false(size(ok));
            end
        end

        denom = sum(ok);

        if denom > 0
            sig = (p(ok) < alpha/2) | (p(ok) > 1 - alpha/2);
            nSig = sum(sig);
            frac = nSig / denom;
        else
            nSig = 0;
            frac = NaN;
        end

        % --- Robust distribution summaries from stats ---
        sn = tn + "_stat";
        if ismember(sn, allVars)
            s = R.(sn);
            s = s(regMask);
            s = s(~isnan(p));
            s = s(ok);

            nValidStat = numel(s);

            if nValidStat > 0
                medStat = median(s);
                iqrVal  = iqr(s);
                madVal  = mad(s, 1);  % median absolute deviation from median
            else
                medStat = NaN;
                iqrVal  = NaN;
                madVal  = NaN;
            end
        else
            nValidStat = 0;
            medStat = NaN;
            iqrVal = NaN;
            madVal = NaN;
        end

        regionCol(rr)      = rID;
        testNameCol(rr)    = tn;

        fracCol(rr)        = frac;
        nValidCol(rr)      = denom;
        nSigCol(rr)        = nSig;
        nROIsCol(rr)       = nROIs;

        medianStatCol(rr)  = medStat;
        iqrStatCol(rr)     = iqrVal;
        madStatCol(rr)     = madVal;
        nValidStatCol(rr)  = nValidStat;
    end
end

Tsess = table(subjCol, dayCol, dRelCol, sessIDCol, fovCol, regionCol, testNameCol, ...
              fracCol, nValidCol, nSigCol, nROIsCol, ...
              medianStatCol, iqrStatCol, madStatCol, nValidStatCol, ...
              'VariableNames', {'subjStr','day','dRelBias','sessionID','fov','regionId','testName', ...
                                'fracSig','nValid','nSig','nROIs', ...
                                'medianStat','iqrStat','madStat','nValidStat'});

end