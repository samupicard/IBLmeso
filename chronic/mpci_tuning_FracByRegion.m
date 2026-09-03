function Tsess = mpci_tuning_FracByRegion(sessionPath, alpha, opt)
% mpci_tuning_FracByRegion
% Compute per-session per-region per-test fraction significant using ALL ROIs
% filtered by mpciROITypes==1 (no chronic filtering).
%
% Output Tsess is LONG format with rows:
%   subjStr, day, dRelBias, sessionID, fov(="ALL"), regionId, testName,
%   fracSig, nValid, nSig, nROIs
%
% Notes:
% - Variable test set across FOVs is supported; missing tests => NaN.
% - Two-tailed criterion: p < alpha/2 OR p > 1-alpha/2.
% - Region filtering via expand_region_ids(opt) applies to regionId.
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
R = table();                 % per-ROI table (regionId + test columns)
allTestNames = string.empty(1,0);

for f = 1:numel(fovs)
    fovPath = fullfile(fovs(f).folder, fovs(f).name);

    pFile    = fullfile(fovPath, 'mpciROIs.taskTunedP.tsv');
    regFile  = fullfile(fovPath, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy');
    typeFile = fullfile(fovPath, 'mpciROIs.mpciROITypes.npy');

    if ~(isfile(pFile) && isfile(regFile) && isfile(typeFile))
        continue
    end

    % Load P-values
    pTab = readtable(pFile, 'FileType','text', 'Delimiter','\t','VariableNamingRule','preserve');
    testNames = string(pTab.Properties.VariableNames);
    P = table2array(pTab);
    n = size(P,1);

    % Union tests
    newTests = setdiff(testNames, allTestNames, 'stable');
    if ~isempty(newTests)
        % add NaN columns to already accumulated R
        for nt = 1:numel(newTests)
            R.(newTests(nt)) = nan(height(R), 1);
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

    % Init all tests in this session union as NaN
    for k = 1:numel(allTestNames)
        T.(allTestNames(k)) = nan(height(T), 1);
    end

    % Fill only tests present in this file
    for k = 1:numel(testNames)
        v = testNames(k);
        tmp = P(:,k);
        T.(v) = tmp(keep);
    end

    % Robust append union columns (align T and R)
    if isempty(R)
        R = T;
    else
        missingInR = setdiff(T.Properties.VariableNames, R.Properties.VariableNames, 'stable');
        for m = 1:numel(missingInR)
            v = missingInR{m};
            R.(v) = nan(height(R), 1);
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

% Determine test columns
testCols = setdiff(string(R.Properties.VariableNames), "regionId", 'stable');

% Compute per-region per-test fractions (long format)
regionList = unique(double(R.regionId));
nRegions = numel(regionList);
nTests = numel(testCols);

rows = nRegions * nTests;

subjCol      = repmat(subjStr, rows, 1);
dayCol       = repmat(day, rows, 1);
dRelCol      = repmat(dRelBias, rows, 1);
sessIDCol    = repmat(sessionID, rows, 1);
fovCol       = repmat("ALL", rows, 1);
regionCol    = nan(rows, 1);
testNameCol  = strings(rows, 1);
fracCol      = nan(rows, 1);
nValidCol    = zeros(rows, 1);
nSigCol      = zeros(rows, 1);
nROIsCol     = zeros(rows, 1);

rr = 0;
for ri = 1:nRegions
    rID = regionList(ri);
    regMask = (double(R.regionId) == rID);
    nROIs = sum(regMask);

    for ti = 1:nTests
        rr = rr + 1;

        p = R.(testCols(ti));
        p = p(regMask);

        ok = ~isnan(p);
        denom = sum(ok);

        if denom > 0
            sig = (p(ok) < alpha/2) | (p(ok) > 1 - alpha/2);
            nSig = sum(sig);
            frac = nSig / denom;
        else
            nSig = 0;
            frac = NaN;
        end

        regionCol(rr)   = rID;
        testNameCol(rr) = testCols(ti);
        fracCol(rr)     = frac;
        nValidCol(rr)   = denom;
        nSigCol(rr)     = nSig;
        nROIsCol(rr)    = nROIs;
    end
end

Tsess = table(subjCol, dayCol, dRelCol, sessIDCol, fovCol, regionCol, testNameCol, ...
              fracCol, nValidCol, nSigCol, nROIsCol, ...
              'VariableNames', {'subjStr','day','dRelBias','sessionID','fov','regionId','testName', ...
                                'fracSig','nValid','nSig','nROIs'});

end