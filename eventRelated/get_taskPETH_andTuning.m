function PETH_struct = get_taskPETH_andTuning(datpath, varargin)
%GET_TASKPETH_ANDTUNING Compute/save task PETHs and task-tuning statistics.
%
% PETH_struct = get_taskPETH_andTuning(datpath, params_all, Name,Value,...)
%
% The same params_all struct array drives both analyses. Fields used by the
% PETH calculation and tuning test are listed below. Missing optional fields
% are assigned sensible defaults.
%
% Common fields:
%   activity_type      default 'deconv'
%   trialTypeField
%   trialTypeVals       [] -> infer from trials table
%   trialTypeFilter     default ''
%   evnt
%   twin_all            full PETH window
%   twin_bl             default 'none'
%   twin_ev             event-response window
%   nTrialsMin          default 20
%   minTrialsPerCond    default 10
%
% PETH-specific:
%   condFields          conditioning fields for getMeanResps_allROIs
%                       ('cndFields' is accepted as a legacy alias)
%   minTrialsPerCombo   default 2
%
% Tuning-specific:
%   stat_to_use         default 'ccu'
%   nComp               default 1
%   nTrialsToKeep       default false
%   statCV              'odd', 'even', or '' (legacy field 'cv' accepted)
%   pthresh             default 0.05
%   nPseudoSessions     default from Name/Value option 'nPseudos'
%
% Name/Value options:
%   'saveToFOV'         false by default
%   'pethCV'            save odd/even PETHs (default false)
%   'cv'                legacy alias for pethCV
%   'nPseudos'          default 199
%   'newPseudos'        default false; preserves legacy custom pseudo mode
%   'overwriteExisting' default false; overwrite existing taskTuned TSV columns
%   'overwriteTuningTSVs'default false; re-write TSV entirely
%   'stat_names'        [] -> infer_stat_names(params_all)
%   'computePETH'       default true
%   'computeTuning'     default true
%
% Outputs written to each alf/FOV_XX folder when saveToFOV=true:
%   mpciROIs.PETHavgNorm_<event>_<trialType>.npy
%   mpciROIs.PETHavgNorm_<event>_<trialType>_odd.npy   (if pethCV)
%   mpciROIs.PETHavgNorm_<event>_<trialType>_even.npy  (if pethCV)
%   mpciROIs.PETHavgNorm_<event>_<trialType>.timeValues.npy
%   mpciROIs.PETHavgNorm_<event>_<trialType>.conditionValues.npy
%   mpciROIs.taskTunedP.tsv
%   mpciROIs.taskTunedStat.tsv
%
% Event extraction is performed once per parameter set with get_evoked_full.
% Its Evk output drives task-tuning statistics, while Dat drives the PETHs.
% A fifth output from get_evoked_full (badTrialsEvk) distinguishes trials
% invalid for the response window from trials invalid for the full PETH.

%% optional positional params_all
if ~isempty(varargin) && isstruct(varargin{1})
    params_all = varargin{1};
    varargin(1) = [];
else
    params_all = defaultParams();
end

%% options
p = inputParser;
p.addParameter('saveToFOV', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('pethCV', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('cv', [], @(x)isempty(x) || islogical(x) || isnumeric(x)); % legacy alias
p.addParameter('nPseudos', 199, @(x)isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('newPseudos', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('overwriteExisting', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('overwriteTuningTSVs', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('stat_names', [], @(x)isempty(x) || iscell(x) || isstring(x));
p.addParameter('computePETH', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('computeTuning', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

saveToFOV         = logical(p.Results.saveToFOV);
pethCV            = logical(p.Results.pethCV);
if ~isempty(p.Results.cv)
    pethCV = logical(p.Results.cv);
end
nPseudosDefault    = p.Results.nPseudos;
newPseudos         = logical(p.Results.newPseudos);
overwriteExisting = logical(p.Results.overwriteExisting);
overwriteTuningTSVs = logical(p.Results.overwriteTuningTSVs);
computePETH        = logical(p.Results.computePETH);
computeTuning      = logical(p.Results.computeTuning);

params_all = normalizeParams(params_all, nPseudosDefault);

if isempty(p.Results.stat_names)
    stat_names = infer_stat_names(params_all);
else
    stat_names = cellstr(p.Results.stat_names);
end
assert(numel(stat_names) == numel(params_all), ...
    'stat_names length must match params_all length.');

%% path/session metadata
splitPath = split(datpath, filesep);
subj = splitPath{end-2};
day  = splitPath{end-1};
sess = splitPath{end};

fprintf('\n%s:\n', datpath);

%% decide which tuning tests actually need computing
runTuningMask = false(1, numel(params_all));
if computeTuning
    runTuningMask(:) = true;

    if saveToFOV && ~overwriteExisting && ~overwriteTuningTSVs
        runTuningMask = findMissingTuningColumns(datpath, stat_names);

        if ~any(runTuningMask)
            fprintf('All task-tuning columns already exist; tuning calculation will be skipped.\n');
        else
            fprintf('Will compute %d/%d tuning tests.\n', ...
                sum(runTuningMask), numel(runTuningMask));
        end
    end
end

%% load trials once
fprintf('Getting event timings... ');
trialsT = IBL_loadTrialsTable(datpath, 'sync', 'timeline');
if isempty(trialsT)
    warning('No extracted trials found. Trying raw Bpod events...');
    try
        trialsT = IBL_generateTrialsTable(datpath);
    catch
        trialsT = array2table([]);
    end
end
if isempty(trialsT)
    fprintf('No trials found, skipping this session.\n');
    PETH_struct = [];
    return
end
trialsT.stimSide = sign(trialsT.contrastDiff);
trialsT.stimSide(trialsT.contrastDiff == 0) = 0;
fprintf('Done!\n');

%% load neural data once
frameQCPath = fullfile(datpath, 'alf', 'FOV_00', 'mpci.mpciFrameQC.npy');
if exist(frameQCPath, 'file')
    frameQC = readNPY(frameQCPath);
    if any(frameQC ~= 0)
        warning('Badframes / non-zero frameQC found. Trials overlapping them will be ignored.');
    end
end

fprintf('Loading 2PI traces...');
Fall = IBL_loadMesoData(subj, day, sess, 'trace', 'spks', 'fast', true);
if isempty(Fall.tr) || isempty(Fall.time)
    warning('Incomplete datasets, skipping this session.');
    PETH_struct = [];
    return
end
if size(Fall.tr,2) ~= numel(Fall.time)
    warning('Unequal frame counts in ROIActivity (%d) and times (%d). Skipping session.', ...
        size(Fall.tr,2), numel(Fall.time));
    PETH_struct = [];
    return
end

sig        = Fall.tr;
frameTimes = Fall.time(:);
fov        = Fall.fov(:);
idx        = Fall.idx(:);
Fs         = 1 / median(diff(frameTimes));
badframes  = find(Fall.frameQC ~= 0);

if isfield(Fall, 'timeshift')
    tshifts = Fall.timeshift;
else
    tshifts = [];
end
nROIs = size(sig,1);

%% sampling guardrail: upsample if minimum time window < imaging timestep
[sig, frameTimes, badframes, Fs] = maybeUpsample( ...
    sig, frameTimes, badframes, Fs, params_all);

%% only load/generate pseudos if at least one tuning test will run
pseudoSets = struct();
if any(runTuningMask)
    maxPseudosNeeded = max([params_all(runTuningMask).nPseudoSessions]);
    pseudoSets = loadOrGeneratePseudos(datpath, trialsT, maxPseudosNeeded);
end

%% initialise output
PETH_struct = struct();
PETH_struct.meta = struct( ...
    'datpath', datpath, ...
    'subject', subj, ...
    'date', day, ...
    'session', sess, ...
    'Fs', Fs, ...
    'nROIs', nROIs, ...
    'nFrames', size(sig,2));
PETH_struct.roi = struct('fov', fov, 'idx', idx);
PETH_struct.params_all = params_all;
PETH_struct.entries = repmat(struct(), 1, numel(params_all));
PETH_struct.out = repmat(struct(), 1, numel(params_all));

statsMat   = nan(nROIs, numel(params_all));
psMat      = nan(nROIs, numel(params_all));
hsMat      = false(nROIs, numel(params_all));
didCompute = false(1, numel(params_all));

% Open one thread pool only if tuning is actually needed.
if any(runTuningMask)
    pobj = gcp('nocreate');
    if isempty(pobj)
        parpool('threads');
    end
end

%% shared parameter loop
for iParam = 1:numel(params_all)
    params = params_all(iParam);
    [evnt, evnt_nm] = resolveEventName(params.evnt);

    if height(trialsT) < params.nTrialsMin
        fprintf('%s/%s: fewer than %d total trials; skipping parameter set.\n', ...
            evnt, params.trialTypeField, params.nTrialsMin);
        continue
    end

    if isempty(params.trialTypeVals)
        if ~ismember(params.trialTypeField, trialsT.Properties.VariableNames)
            error('Trial table has no variable "%s".', params.trialTypeField);
        end
        x = trialsT.(params.trialTypeField);
        trialTypeVals = unique(x(~isnan(x)))';
    else
        trialTypeVals = params.trialTypeVals;
    end

    % Extract once. Evk is used for tuning; Dat is used for mean PETHs.
    % badTrialsPETH reflects validity of the full twin_all window, whereas
    % badTrialsStat reflects validity of twin_ev/twin_bl only.
    allResps = [];
    allPETH_raw = [];
    T = [];
    badTrialsPETH = [];
    badTrialsStat = [];

    if computePETH || runTuningMask(iParam)
        fprintf('Extracting %s PETHs time-locked to %s... ', params.trialTypeField, evnt);
        try
            [allResps, allPETH_raw, T, badTrialsPETH, badTrialsStat] = get_evoked_full( ...
                sig, frameTimes, badframes, trialsT, evnt_nm, ...
                params.twin_ev, params.twin_bl, params.twin_all, tshifts);
        catch ME
            warning('Event extraction failed for %s/%s:\n%s', ...
                evnt, params.trialTypeField, getReport(ME,'basic'));
        end
    end

    %% PETH branch
    if computePETH && ~isempty(allPETH_raw)
        entry = makePETHEntry(allPETH_raw, T, badTrialsPETH, sig, frameTimes, ...
            trialsT, params, trialTypeVals, pethCV);
        entry.activity_type = params.activity_type;
        entry.evnt = evnt;
        entry.evnt_nm = evnt_nm;
        entryFields = fieldnames(entry);
        for iField = 1:numel(entryFields)
            fn = entryFields{iField};
            PETH_struct.entries(iParam).(fn) = entry.(fn);
        end
    else
        PETH_struct.entries(iParam).activity_type = params.activity_type;
        PETH_struct.entries(iParam).evnt = evnt;
        PETH_struct.entries(iParam).evnt_nm = evnt_nm;
    end

    %% tuning branch
    if runTuningMask(iParam) && ~isempty(allResps)
        if newPseudos
            fprintf('Generating custom pseudos for %s... ', stat_names{iParam});
            pseudoSessions = repmat(struct('probabilityLeft',[], 'contrastDiff',[], ...
                'choice',[], 'feedbackType',[]), 1, params.nPseudoSessions);
            for iP = 1:params.nPseudoSessions
                pseudoSessions(iP) = IBL_genSession_feedback0(trialsT);
            end
            fprintf('Done!\n');
        else
            pseudoSessions = selectPseudoSet(pseudoSets, params.trialTypeField);
            if isempty(pseudoSessions)
                warning('No suitable pseudosessions for %s; skipping tuning test.', stat_names{iParam});
                continue
            end
        end

        fprintf('Tuning test for %s %.1f v. %.1f: ', params.trialTypeField, params.compareVals(1), params.compareVals(2));
        [statReal, statP, statH, statsAll, ok] = computeTuningTest( ...
            allResps, badTrialsStat, trialsT, pseudoSessions, params, nROIs);

        if ~ok
            warning('Tuning test %s failed its real-session validity criterion.', stat_names{iParam});
            continue
        end
        if all(isnan(statReal)) || all(isnan(statP))
            warning('Tuning test %s produced all-NaN output.', stat_names{iParam});
            continue
        end

        PETH_struct.out(iParam).stats = statsAll;
        PETH_struct.out(iParam).ps = statP;
        PETH_struct.out(iParam).hs = statH;

        statsMat(:,iParam) = statReal;
        psMat(:,iParam) = statP;
        hsMat(:,iParam) = statH;
        didCompute(iParam) = true;
    end
end

%% save all outputs in one FOV pass
if saveToFOV
    fprintf('Saving PETHs and tuning results... ');
    saveAllToFOV(datpath, fov, PETH_struct.entries, pethCV, ...
        statsMat, psMat, didCompute, stat_names, overwriteExisting, overwriteTuningTSVs);
    fprintf('Done!\n');
end

end


%% ------------------------------------------------------------------------

function runMask = findMissingTuningColumns(datpath, stat_names)
% Run a test if its column is absent from at least one FOV's taskTunedP.tsv.
runMask = false(1, numel(stat_names));
fovDirs = dir(fullfile(datpath, 'alf', 'FOV_*'));

if isempty(fovDirs)
    runMask(:) = true;
    return
end

for d = 1:numel(fovDirs)
    if ~fovDirs(d).isdir, continue; end
    tsvPath = fullfile(fovDirs(d).folder, fovDirs(d).name, 'mpciROIs.taskTunedP.tsv');
    if ~exist(tsvPath, 'file')
        runMask(:) = true;
        continue
    end
    try
        Tprev = readtable(tsvPath, 'FileType','text', 'Delimiter','\t', ...
            'PreserveVariableNames', true);
        runMask = runMask | ~ismember(stat_names, Tprev.Properties.VariableNames);
    catch
        runMask(:) = true;
    end
end
end


function [sig, frameTimes, badframes, Fs] = maybeUpsample(sig, frameTimes, badframes, Fs, params_all)
allBinSizes = nan(1, numel(params_all));
for i = 1:numel(params_all)
    x = params_all(i).twin_ev;
    if isnumeric(x) && numel(x) == 2
        allBinSizes(i) = abs(diff(x));
    elseif iscell(x) && ~isempty(x)
        allBinSizes(i) = min(cellfun(@(y)abs(diff(y)), x));
    end
end

if all(isnan(allBinSizes)), return; end
minBin = min(allBinSizes(~isnan(allBinSizes)));
dt = 1/Fs;
if dt <= minBin, return; end

newDt = floor(20*(0.99*dt))/20;
if newDt <= 0, return; end
warning(['Sampling interval (%.4fs) > min twin_ev bin (%.4fs). ' ...
    'Resampling to %.4fs without interpolating across gaps.'], dt, minBin, newDt);

gapIdx = find(diff(frameTimes) > 3*dt);
segStarts = [1; gapIdx(:)+1];
segEnds = [gapIdx(:); numel(frameTimes)];
newT = cell(numel(segStarts),1);
newS = cell(numel(segStarts),1);
for i = 1:numel(segStarts)
    ix = segStarts(i):segEnds(i);
    if numel(ix) < 2, continue; end
    tSeg = frameTimes(ix);
    tNew = (tSeg(1):newDt:tSeg(end))';
    if isempty(tNew), continue; end
    newT{i} = tNew;
    newS{i} = interp1(tSeg, sig(:,ix)', tNew, 'linear')';
end
keep = ~cellfun(@isempty,newT);
if ~any(keep), return; end

oldTimes = frameTimes;
oldBad = badframes;
frameTimes = vertcat(newT{keep});
sig = cat(2,newS{keep});
Fs = 1/median(diff(frameTimes));

badframes = [];
if ~isempty(oldBad)
    oldBadTimes = oldTimes(oldBad);
    [nearestIdx, nearestDist] = knnsearch(frameTimes, oldBadTimes);
    badframes = unique(nearestIdx(nearestDist <= newDt/2))';
end
for i = 1:numel(gapIdx)
    a = oldTimes(gapIdx(i));
    b = oldTimes(gapIdx(i)+1);
    badframes = unique([badframes, find(frameTimes >= a-newDt & frameTimes <= b+newDt)']); %#ok<AGROW>
end
end



function [evnt, evnt_nm] = resolveEventName(evntIn)
evntIn = char(evntIn);
if endsWith(evntIn, '_times')
    evnt_nm = evntIn;
    evnt = extractBefore(evntIn, strlength(evntIn)-strlength('_times')+1);
    evnt = char(evnt);
else
    evnt = evntIn;
    evnt_nm = [evntIn '_times'];
end
end



function entry = makePETHEntry(allPETH_raw, T, badTrials, sig, frameTimes, trialsT, params, trialTypeVals, pethCV)
% Normalize, filter, recode contrastDiff to side, then average by condition.
taskTimes = [trialsT.intervals_0(1), trialsT.intervals_1(end)];
taskFrames = frameTimes > taskTimes(1) & frameTimes < taskTimes(2);
prctiles = prctile(sig(:,taskFrames), [20,99], 2);
allNorm = normalize_minmax(allPETH_raw, prctiles(:,1), prctiles(:,2));

trialKeep = ~badTrials(:);
filterMask = evalTrialTypeFilter(trialsT, params.trialTypeFilter);
trialKeep = trialKeep & filterMask(:);
trialsQ = trialsT(trialKeep,:);
allNorm = allNorm(:,trialKeep,:);

entry = struct();
entry.trialTypeField = params.trialTypeField;
entry.trialTypeVals = trialTypeVals;
entry.trialTypeFilter = params.trialTypeFilter;
entry.compareVals    = params.compareVals;
entry.twin_all = params.twin_all;
entry.twin_bl = params.twin_bl;
entry.twin_ev = params.twin_ev;
entry.T = T(:)';
entry.badTrials = badTrials(:);
entry.minTrialsPerCond = params.minTrialsPerCond;

if isempty(trialsQ)
    warning('No trials remain after QC/filter for PETH %s.', params.trialTypeField);
    return
end

x = trialsQ.(params.trialTypeField);
nTrialsPerCond = arrayfun(@(v)sum(x==v), trialTypeVals);
insufficient = nTrialsPerCond < params.minTrialsPerCond;
entry.nTrialsPerCond = nTrialsPerCond;

PETH = getMeanResps_allROIs(allNorm, trialsQ, params.trialTypeField, ...
    trialTypeVals, params.condFields, params.minTrialsPerCombo);
PETH(:,insufficient,:) = NaN;
entry.PETH = PETH;

if pethCV
    odd = 1:2:height(trialsQ);
    even = 2:2:height(trialsQ);
    entry.PETH_odd = getMeanResps_allROIs(allNorm(:,odd,:), trialsQ(odd,:), ...
        params.trialTypeField, trialTypeVals, params.condFields, params.minTrialsPerCombo);
    entry.PETH_even = getMeanResps_allROIs(allNorm(:,even,:), trialsQ(even,:), ...
        params.trialTypeField, trialTypeVals, params.condFields, params.minTrialsPerCombo);
    entry.PETH_odd(:,insufficient,:) = NaN;
    entry.PETH_even(:,insufficient,:) = NaN;
end
end


function [statReal, statP, statH, stats, realOK] = computeTuningTest(allResps, badTrials, trialsT, pseudoSessions, params, nROIs)
nTrials = size(allResps,2);
baseValid = true(1,nTrials);
baseValid(max(1,nTrials-19):nTrials) = false;
baseValid(any(isnan(allResps),1)) = false;
baseValid(all(allResps==0,1)) = false;
baseValid(badTrials(:)') = false;

if strcmpi(params.statCV,'odd')
    baseValid(2:2:end) = false;
elseif strcmpi(params.statCV,'even')
    baseValid(1:2:end) = false;
end

nPseudoUse = min(params.nPseudoSessions, numel(pseudoSessions));
nSess = 1+nPseudoUse;
sessionsT = cell(1,nSess);
sessionsT{1} = trialsT;
for k = 1:nPseudoUse
    sessionsT{k+1} = struct2table(pseudoSessions(k));
end

stats = nan(nROIs,nSess);
ok = false(1,nSess);
err = strings(1,nSess);

trialTypeField = params.trialTypeField;
trialTypeFilter = params.trialTypeFilter;
condFields = params.condFields;
nTrialsMin = params.nTrialsMin;
nTrialsToKeep = params.nTrialsToKeep;
stat_to_use = params.stat_to_use;

fprintf('Parallel-computing %s over %d real/pseudo sessions... ', stat_to_use, nSess);
parfor s = 1:nSess
    try
        sessT = sessionsT{s};
        sessT.stimSide = sign(sessT.contrastDiff);
        sessT.stimSide(sessT.contrastDiff == 0) = 0;

        compareVals = params.compareVals;
        
        if numel(compareVals) ~= 2
            error('compareVals must contain exactly two values.');
        end
        
        trialVals = sessT.(trialTypeField);
        
        group1 = trialVals == compareVals(1);
        group2 = trialVals == compareVals(2);
        
        filterMask = evalTrialTypeFilter(sessT, trialTypeFilter);
        
        group1 = group1(:)' & filterMask(:)';
        group2 = group2(:)' & filterMask(:)';
        
        compMask = group1 | group2;
        validMask = baseValid & compMask;

        if nTrialsMin && sum(validMask) < nTrialsMin
            continue
        end
        if nTrialsToKeep && sum(validMask) > nTrialsToKeep
            tmp = false(size(validMask));
            ii = find(validMask, nTrialsToKeep, 'first');
            tmp(ii) = true;
            validMask = tmp;
        end

        iC1 = group1(validMask);
        iC2 = group2(validMask);
        validResps = allResps(:,validMask);

        switch stat_to_use
            case 'dprime'
                test_stat = getdprimes(validResps,[iC2;iC1]);
            case 'dprime_cv'
                test_stat = getdprimes(validResps,[iC2;iC1],5);
            case 'auc'
                test_stat = getAUCs(validResps,[iC2;iC1]);
            case 'auc_cv'
                test_stat = getAUCs(validResps,[iC2;iC1],5);
            case 'meandiff'
                test_stat = getmeandiffs(validResps,[iC2;iC1]);
            case 'ccMean'
                test_stat = getccMeans(validResps,sessT(validMask,:),trialTypeField,condFields);
            case 'ccMeanDiff'
                test_stat = getccMeanDiffs(validResps,sessT(validMask,:),iC1,trialTypeField,condFields);
            case 'ccu'
                test_stat = getccUs(validResps,sessT(validMask,:),iC1,trialTypeField,condFields);
            case 'ccMI'
                test_stat = getccMIs(validResps,sessT(validMask,:),iC1,trialTypeField,condFields);
            otherwise
                test_stat = nan(nROIs,1);
        end
        stats(:,s) = test_stat(:);
        ok(s) = true;
    catch ME
        err(s) = string(ME.message);
    end
end
fprintf('Done!\n');

bad = find(err ~= "");
if ~isempty(bad)
    disp([bad(:), err(bad).']);
end

realOK = ok(1);
statReal = stats(:,1);
statP = nan(nROIs,1);
statH = false(nROIs,1);
if ~realOK, return; end

goodPseudo = find(ok(2:end));
if isempty(goodPseudo)
    warning('No valid pseudosessions passed the trial criteria.');
    realOK = false;
    return
end
stats = [stats(:,1), stats(:,1+goodPseudo)];
for i = 1:nROIs
    statP(i) = 0.01 * comp_percentile(stats(i,2:end), stats(i,1));
end
statH = statP > 1-params.pthresh/2 | statP < params.pthresh/2;
end


function saveAllToFOV(datpath, fov, entries, pethCV, ...
    statsMat, psMat, didCompute, stat_names, ...
    overwriteExisting, overwriteTuningTSVs)

uFov = unique(fov(:))';
makeSafe = @(s)regexprep(char(s),'[^A-Za-z0-9_\-]','');

for k = 1:numel(uFov)

    thisFov = uFov(k);
    fovFolder = fullfile(datpath,'alf',sprintf('FOV_%02d',thisFov));

    if ~exist(fovFolder,'dir')
        warning('FOV folder not found: %s', fovFolder);
        continue
    end

    roiMask = fov == thisFov;
    if ~any(roiMask)
        continue
    end

    %% PETH files
    for i = 1:numel(entries)

        e = entries(i);

        if ~isfield(e,'PETH') || isempty(e.PETH)
            continue
        end

        fName = sprintf( ...
            'PETHavgNorm_%s_%s', ...
            makeSafe(e.evnt), ...
            makeSafe(e.trialTypeField));

        baseROIs = fullfile(fovFolder, ['mpciROIs.' fName]);
        base = fullfile(fovFolder, fName);

        writeNPY(e.PETH(roiMask,:,:), [baseROIs '.npy']);

        if pethCV && isfield(e,'PETH_odd') && isfield(e,'PETH_even')
            writeNPY(e.PETH_odd(roiMask,:,:),  [baseROIs '_odd.npy']);
            writeNPY(e.PETH_even(roiMask,:,:), [baseROIs '_even.npy']);
        end

        if isfield(e,'T') && ~isempty(e.T)
            writeNPY(e.T(:), [base '.timeValues.npy']);
        end

        if isfield(e,'trialTypeVals') && ~isempty(e.trialTypeVals)
            writeNPY(e.trialTypeVals(:), [base '.conditionValues.npy']);
        end
    end

    %% taskTuned TSVs
    successNames = stat_names(didCompute);

    if isempty(successNames)
        continue
    end

    TnewP = array2table( ...
        psMat(roiMask,didCompute), ...
        'VariableNames',successNames);

    TnewS = array2table( ...
        statsMat(roiMask,didCompute), ...
        'VariableNames',successNames);

    outP = fullfile(fovFolder,'mpciROIs.taskTunedP.tsv');
    outS = fullfile(fovFolder,'mpciROIs.taskTunedStat.tsv');

    if overwriteTuningTSVs

        % Completely replace existing TSVs with the results from this run.
        writetable(TnewP,outP,'FileType','text','Delimiter','\t');
        writetable(TnewS,outS,'FileType','text','Delimiter','\t');

    else

        mergeTuningTables( ...
            outP, outS, ...
            TnewP, TnewS, ...
            overwriteExisting);

    end
end
end

function mergeTuningTables(outP,outS,TnewP,TnewS,overwriteExisting)
if ~exist(outP,'file') || ~exist(outS,'file')
    writetable(TnewP,outP,'FileType','text','Delimiter','\t');
    writetable(TnewS,outS,'FileType','text','Delimiter','\t');
    return
end

try
    ToldP = readtable(outP,'FileType','text','Delimiter','\t','PreserveVariableNames',true);
    ToldS = readtable(outS,'FileType','text','Delimiter','\t','PreserveVariableNames',true);
catch ME
    warning('Could not read existing tuning TSVs (%s); writing new results only.',ME.message);
    writetable(TnewP,outP,'FileType','text','Delimiter','\t');
    writetable(TnewS,outS,'FileType','text','Delimiter','\t');
    return
end

if height(ToldP) ~= height(TnewP) || height(ToldS) ~= height(TnewS) || ...
        ~isequal(ToldP.Properties.VariableNames,ToldS.Properties.VariableNames)
    warning('Existing tuning TSVs are incompatible with current ROI rows; writing new results only.');
    writetable(TnewP,outP,'FileType','text','Delimiter','\t');
    writetable(TnewS,outS,'FileType','text','Delimiter','\t');
    return
end

oldNames = ToldP.Properties.VariableNames;
newNames = TnewP.Properties.VariableNames;
ToutP = ToldP;
ToutS = ToldS;

if overwriteExisting
    common = intersect(oldNames,newNames,'stable');
    for i = 1:numel(common)
        nm = common{i};
        ToutP.(nm) = TnewP.(nm);
        ToutS.(nm) = TnewS.(nm);
    end
end

add = setdiff(newNames,oldNames,'stable');
if ~isempty(add)
    ToutP = [ToutP,TnewP(:,add)];
    ToutS = [ToutS,TnewS(:,add)];
end

writetable(ToutP,outP,'FileType','text','Delimiter','\t');
writetable(ToutS,outS,'FileType','text','Delimiter','\t');
end
