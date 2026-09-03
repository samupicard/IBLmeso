function PETH_struct = get_taskTunedROIs(datpath, varargin)
%
% get_taskTunedROIs
%
% loads trial data and all ROI traces from all FOVs in a IBL meso session
% computes time-averaged PETHs relative to some pre-defined events
% for real session and a set of pseudo-sessions, computes test statistic for
% pre-defined set of task variable comparisons
% saves a matrix of test-statistics and a matrix of p-values for each ROI 
% and each comparison
%
% Samuel Picard (Feb 2026)


%% options

% check for positional params_all struct
if ~isempty(varargin) && isstruct(varargin{1})
    params_all = varargin{1};
    varargin(1) = []; % remove it so inputParser doesn't see it
else
    params_all = [];
end

p = inputParser;
p.addParameter('saveToFOV', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('nPseudos', 199, @(x)isnumeric(x));
p.addParameter('newPseudos', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('overwriteExisting', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('stat_names', [], @(x) isempty(x) || iscell(x) || isstring(x));
p.parse(varargin{:});

saveToFOV         = logical(p.Results.saveToFOV);
nPseudoSessions   = p.Results.nPseudos;
newPseudos        = p.Results.newPseudos;
overwriteExisting = logical(p.Results.overwriteExisting);

%define the event types, windows etc to use for the PETHs
if isempty(params_all)
    params_all = struct( ...
        'activity_type',   'deconv', ...
        'trialTypeField',  {'contrastDiff','contrastDiff','choice','feedbackType','probabilityLeft'}, ...
        'trialTypeVals',   {[-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1],[-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1], [-1,1], [-1,1], [0.2,0.5,0.8]}, ...
        'trialTypeFilter', {'', '', '', '', ''}, ...
        'condFields',      {''}, ...
        'evnt',            {'stimOn','stimOn','choiceMovement','feedback','stimOn'}, ...
        'twin_all',        {[-1,3],[-1,3],[-2.5,1.5],[-1.5,2.5],[-2,2]}, ...
        'twin_bl',         'none', ...
        'twin_ev',         {[0,0.4],[0,0.4],[-0.4,0],[0,0.4],[-0.5,-0.1]}, ...
        'stat_to_use',     {'ccu'}, ...
        'nComp',           {1,5,1,1,1}, ...
        'nTrialsMin',      20, ...
        'nTrialsToKeep',   false, ...
        'minTrialsPerCond',10, ...
        'cv',              'even', ...
        'pthresh',         0.05, ...
        'nPseudoSessions', 199 ...
        );
end

% params_all = struct( ...
%     'activity_type',   'deconv', ...
%     'trialTypeField',  {'contrastDiff','contrastDiff','choice','feedbackType'}, ...
%     'trialTypeVals',   {[-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1],[-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1], [-1,1], [-1,1]}, ...
%     'trialTypeFilter', {'', '', '', ''}, ...
%     'condFields',      {''}, ...
%     'evnt',            {'stimOn','stimOn','choiceMovement','feedback'}, ...
%     'twin_all',        {[-1,3],[-1,3],[-2.5,1.5],[-1.5,2.5]}, ...
%     'twin_bl',         'none', ...
%     'twin_ev',         {[0,0.15],[0,0.15],[-0.15,0],[0,0.15]}, ...
%     'stat_to_use',     {'ccu'}, ...
%     'nComp',           {1,5,1,1}, ...
%     'nTrialsMin',      20, ...
%     'nTrialsToKeep',   false, ...
%     'minTrialsPerCond',10, ...
%     'cv',              'even', ...
%     'pthresh',         0.05, ...
%     'nPseudoSessions', 199 ...
%     );

% params_all = struct( ...
%     'activity_type',   'deconv', ...
%     'trialTypeField',  {'contrastDiff'}, ...
%     'trialTypeVals',   {[-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1]}, ...
%     'trialTypeFilter', {''}, ...
%     'condFields',      {{'choice'}}, ...
%     'evnt',            {'stimOn'}, ...
%     'twin_all',        {[-1,3]}, ...
%     'twin_bl',         'none', ...
%     'twin_ev',         {[0,0.15]}, ...
%     'stat_to_use',     {'ccMI'}, ...
%     'nComp',           {1}, ...
%     'nTrialsMin',      20, ...
%     'nTrialsToKeep',   false, ...
%     'minTrialsPerCond',10, ...
%     'pthresh',         0.05, ...
%     'nPseudoSessions', 199 ...
%     );

% stat_names: infer by default from params_all
if isempty(p.Results.stat_names)
    stat_names = infer_stat_names(params_all);
else
    stat_names = cellstr(p.Results.stat_names);
end

% sanity check
assert(numel(stat_names) == numel(params_all), 'stat_names length must match params_all length.')

%path logic
splitPath = split(datpath, filesep);
subj = splitPath{end-2};
day = splitPath{end-1};
sess = splitPath{end};

fprintf('\n%s:\n',datpath);

%% determine which tests need to run (based on existing saved TSVs)
runMask = true(1, numel(stat_names));  % default: run everything

if saveToFOV
    if overwriteExisting
        % recompute everything (we'll only overwrite matching columns later)
        runMask = true(1, numel(stat_names));
        fprintf('overwriteExisting=true: will recompute all %d tests.\n', numel(stat_names));
    else
        % Look for existing results in any FOV. We run only missing columns (union across FOVs).
        fovBase = fullfile(datpath, 'alf');
        fovDirs = dir(fullfile(fovBase, 'FOV_*'));
        missingAny = false(1, numel(stat_names));

        for d = 1:numel(fovDirs)
            if ~fovDirs(d).isdir, continue; end

            fovFolder = fullfile(fovBase, fovDirs(d).name);
            tsvPath   = fullfile(fovFolder, 'mpciROIs.taskTunedP.tsv');

            if exist(tsvPath, 'file')
                try
                    Tprev = readtable(tsvPath, 'FileType','text', 'Delimiter','\t', 'PreserveVariableNames',true);
                    prevNames = Tprev.Properties.VariableNames;
                catch
                    prevNames = {};
                end
                missingHere = ~ismember(stat_names, prevNames);
                missingAny  = missingAny | missingHere;
            else
                missingAny = true(1, numel(stat_names));
            end
        end

        runMask = missingAny;

        if ~any(runMask)
            fprintf('All stat_names already present in existing mpciROIs.taskTunedP.tsv files. Nothing new to compute, skipping this session.\n');
            return
        else
            fprintf('Will compute %d/%d tests (only missing columns across FOVs).\n', sum(runMask), numel(runMask));
        end
    end
end

%% load task data

fprintf('Getting event timings... ');
trialsT = IBL_loadTrialsTable(datpath,'sync','timeline');
if isempty(trialsT)
    warning('no extracted trials found. Will try to generate trials table from raw bpod events...');
    %try generating the TrialsTable from the bpod time events
    try
        trialsT = IBL_generateTrialsTable(datpath);
    catch
        trialsT = array2table([]);
    end
end
if isempty(trialsT)
    fprintf('no extracted trials found, skipping this session!\n');
    PETH_struct = [];
    return
    %elseif size(trialsT,1)<params_all(1).nTrialsMin
    %    fprintf('fewer than %d total trials, skipping this session!\n', params_all(1).nTrialsMin);
    %    PETH_struct = [];
    %    return
end
fprintf('Done!\n');

%% load neural data

%load frame QC, check if there are any non-zero, and throw warning if so
frameQC = readNPY(fullfile(datpath,'alf','FOV_00','mpci.mpciFrameQC.npy'));
if any(frameQC~=0)
    warning('Badframes / non-zero frameQC found. We will ignore trials in those periods.\n');
end

%load data from all FOVs into a single struct (this takes ~30-60s)
fprintf('Loading 2PI traces...');
Fall = IBL_loadMesoData(subj,day,sess,'trace','spks','fast',true);

if isempty(Fall.tr) || isempty(Fall.time)
    warning('Incomplete datasets, skipping this session!\n');
    PETH_struct = [];
    return
end

if size(Fall.tr,2) ~= size(Fall.time,1)
    warning(sprintf('Unequal nr of frames in mpci.ROIActivity (%d) and in mpci.times (%d). Skipping this session!\n',size(Fall.tr,2),size(Fall.time,1)));
    PETH_struct = [];
    return
end


sig = Fall.tr;              % deconvolved activity
frameTimes = Fall.time;     % relies on existing mpci.times.npy
Fs = 1/median(diff(Fall.time)); % calculate frameRate from time vector
fov = Fall.fov;             % per-ROI fov label (typically 0..N-1)
idx = Fall.idx;             % per-ROI local index within FOV (if provided)
tshifts = Fall.timeshift;

nROIs = size(sig,1);

%additionally derive the badframes from frameQC
badframes = find(Fall.frameQC~=0); % for now, assume all nonzero frameQC is bad

%% Sampling rate guardrail (upsampling if needed)

% median sampling interval
dt = 1 / Fs;

% collect all twin_ev bin sizes
allBinSizes = nan(1, numel(params_all));
for iP = 1:numel(params_all)
    twin_ev = params_all(iP).twin_ev;
    if isnumeric(twin_ev) && numel(twin_ev) == 2
        allBinSizes(iP) = abs(diff(twin_ev));
    elseif iscell(twin_ev)
        tmp = cellfun(@(x) abs(diff(x)), twin_ev);
        allBinSizes(iP) = min(tmp);
    end
end

%upsample if minimum time window < imaging timestep
minBin = min(allBinSizes(~isnan(allBinSizes)));
if dt > minBin

    % choose minimal required timestep (slightly smaller than bin)
    newDt = floor(20*(0.99*dt))/20; %e.g. 0.15, 0.2, 0.25, etc.

    warning(['Sampling interval (%.4fs) > min twin_ev bin (%.4fs).\n' ...
        'Resampling to dt = %.4fs (Fs = %.2f Hz), without interpolating across gaps.'], ...
        dt, minBin, newDt, 1/newDt);

    % Detect gaps in original acquisition
    gapThresh = 3 * dt;
    gapIdx = find(diff(frameTimes) > gapThresh);

    segStarts = [1; gapIdx(:) + 1];
    segEnds   = [gapIdx(:); numel(frameTimes)];

    newFrameTimesCell = cell(numel(segStarts), 1);
    newSigCell = cell(numel(segStarts), 1);

    for iSeg = 1:numel(segStarts)

        ix = segStarts(iSeg):segEnds(iSeg);

        % Need at least two samples for interpolation
        if numel(ix) < 2
            continue
        end

        tSeg = frameTimes(ix);
        sigSeg = sig(:, ix);

        tNewSeg = (tSeg(1):newDt:tSeg(end))';

        % Avoid creating duplicate or empty segments
        if numel(tNewSeg) < 1
            continue
        end

        newFrameTimesCell{iSeg} = tNewSeg;
        newSigCell{iSeg} = interp1(tSeg(:), sigSeg', tNewSeg, 'linear')';
    end

    keepSeg = ~cellfun(@isempty, newFrameTimesCell);

    frameTimesOld = frameTimes;
    badframesOld = badframes;

    frameTimes = vertcat(newFrameTimesCell{keepSeg});
    sig = cat(2, newSigCell{keepSeg});

    Fs = 1 / median(diff(frameTimes));

    % Remap old bad frames to nearest new frames, but only within tolerance
    badframes = [];

    if ~isempty(badframesOld)
        oldBadTimes = frameTimesOld(badframesOld);

        [nearestIdx, nearestDist] = knnsearch(frameTimes(:), oldBadTimes(:));

        tol = newDt / 2;
        nearestIdx = nearestIdx(nearestDist <= tol);

        badframes = unique(nearestIdx(:))';
    end

    % mark frames adjacent to original gaps as bad, so windows spanning gaps fail
    gapBuffer = newDt;

    for iGap = 1:numel(gapIdx)
        gapStartT = frameTimesOld(gapIdx(iGap));
        gapEndT   = frameTimesOld(gapIdx(iGap) + 1);

        nearGap = frameTimes >= gapStartT - gapBuffer & ...
            frameTimes <= gapEndT   + gapBuffer;

        badframes = unique([badframes, find(nearGap)']);
    end

    fprintf('Resampled from %d to %d frames across %d continuous imaging segments.\n', ...
        numel(frameTimesOld), numel(frameTimes), sum(keepSeg));

end

%% load (or generate) pseudo-sessions

%TO DO properly deal with pseudoSessions for feedback & choice
nPseudos = 300;
pseudoS_flag = false;
fileInfo = dir(fullfile(datpath,'alf','pseudoSessions.mat'));
if ~isempty(fileInfo) && datetime(fileInfo.datenum, 'ConvertFrom', 'datenum') > datetime(2026,06,23)%15)
    fprintf('Loading pseudosessions..');
    load(fullfile(datpath,'alf','pseudoSessions.mat'));
    %check that pseudoSessions and real session have the same nr. of trials
    if size(pseudoSessions_contrast(1).probabilityLeft,1)~=size(trialsT,1)
        pseudoS_flag = true;
        warning('Pre-loaded pseudosessions do not have the same nr of trials as real session! Will regenerate and overwrite.')
    end
    if length(pseudoSessions_contrast)<=250 %strcmp(datpath,'Y:\Subjects\SP054\2024-02-21\001')
        pseudoS_flag = true;
        warning('Not enough pseudosessions available - will regenerate and overwrite.')
    end
else
    pseudoS_flag = true;
end
if pseudoS_flag
    fprintf('Generating pseudosessions..');
    pseudoSessions_bias = struct('probabilityLeft',[],'contrastDiff',[],'choice',[],'feedbackType',[]);
    pseudoSessions_contrast = struct('probabilityLeft',[],'contrastDiff',[],'choice',[],'feedbackType',[]);
    pseudoSessions_choice = struct('probabilityLeft',[],'contrastDiff',[],'choice',[],'feedbackType',[]);
    pseudoSessions_feedback = struct('probabilityLeft',[],'contrastDiff',[],'choice',[],'feedbackType',[]);

    if numel(unique(trialsT.probabilityLeft))==3
        protocolType='bCW';
    else
        protocolType='tCW';
    end

    for iP=1:nPseudos
        if strcmp(protocolType,'bCW')
            pseudoSessions_bias(iP) = IBL_genSession(size(trialsT,1),trialsT); %re-draw blocks and re-draw contrasts
        end
        pseudoSessions_contrast(iP) = IBL_permSession(trialsT,'contrastDiff','pairs'); %randomly permute contrasts (within choice and block)
        pseudoSessions_choice(iP) = IBL_permSession(trialsT,'choice','pairs'); %swap choice in pairs (within contrastDiff and probabilityLeft)
        pseudoSessions_feedback(iP) = IBL_permSession(trialsT,'feedbackType','pairs'); %swap feedback in pairs (within choice and probabilityLeft)
    end
    if strcmp(protocolType,'bCW')
        save(fullfile(datpath,'alf','pseudoSessions.mat'),'pseudoSessions_bias', 'pseudoSessions_contrast', 'pseudoSessions_choice', 'pseudoSessions_feedback');
    else
        save(fullfile(datpath,'alf','pseudoSessions.mat'), 'pseudoSessions_contrast', 'pseudoSessions_choice', 'pseudoSessions_feedback');
    end
end

fprintf('. Done!\n')

%% compute mean PETHs

%initialize PETH struct (session-level metadata + one entry per param)
PETH_struct = struct();
PETH_struct.meta = struct( ...
    'datpath', datpath, ...
    'subject', subj, ...
    'date',    day, ...
    'session', sess, ...
    'Fs',      Fs, ...
    'nROIs',   size(sig,1), ...
    'nFrames', size(sig,2) ...
    );

PETH_struct.roi = struct( ...
    'fov', fov(:), ...
    'idx', idx(:) ...
    );

PETH_struct.params_all = params_all; % keep a copy for provenance
PETH_struct.entries = repmat(struct(), 1, numel(params_all));
PETH_struct.out = repmat(struct(), 1, numel(params_all));

% preallocate outputs for all tests; fill only those that are computed
statsMat = nan(nROIs, numel(params_all));
psMat = nan(nROIs, numel(params_all));
hsMat = false(nROIs, numel(params_all));
didCompute = false(1, numel(params_all));

% Create pool once (optional guard)
pobj = gcp('nocreate');
if isempty(pobj)
    parpool('threads');
end

for iParam = 1:length(params_all)

    if ~runMask(iParam)
        continue
    end

    params = params_all(iParam);

    % append event name with _times (robustly handle if already provided)
    if endsWith(params.evnt, '_times')
        evnt_nm = params.evnt;
        evnt = extractBefore(evnt_nm, strlength(evnt_nm) - strlength('_times') + 1);
        evnt = char(evnt);
    else
        evnt_nm = [params.evnt '_times'];
        evnt = params.evnt;
    end

    % get trialTypeVals
    if isempty(params.trialTypeVals) && isfield(trialsT,params.trialTypeField)
        vals = trialsT.(params.trialTypeField);
        trialTypeVals = unique(vals(~isnan(vals)));
    else
        trialTypeVals = params.trialTypeVals;
    end

    %OPTIONAL: skip if we don't have enough total trials
    if height(trialsT) < params.nTrialsMin
        sprintf('Less than %d trials in total, skipping.\n',params.nTrialsMin);
        continue
    end

    % compute PETHs (raw and time-averaged event triggered responses)
    fprintf('Computing mean responses [%.2f,%.2f]s relative to %s... ', params.twin_ev(1), params.twin_ev(2), evnt);
    try
        [allResps, badTrials] = get_evoked( ...
            sig, frameTimes, badframes, trialsT, evnt_nm, ...
            params.twin_ev, params.twin_bl, tshifts);
    catch ME
        warning('Error computing responses, skipping this type!\n%s\n', getReport(ME, 'basic'));
        %PETH_struct = [];
        continue;
    end

    %CHOOSE WHICH TYPES OF PSEUDOSESSIONS TO USE
    if newPseudos
        pseudoSessions = struct('probabilityLeft',[],'contrastDiff',[],'choice',[],'feedbackType',[]);
        fprintf('Doing custom randomization test (ignoring pre-generated pseudosessions)..');
        for iP=1:nPseudos
            pseudoSessions(iP) = IBL_genSession_feedback0(trialsT); %CHOOSE HERE WHAT TEST TO DO
        end
        fprintf('. Done!\n')
    else
        if strncmpi(params.trialTypeField,'contrast',5)
            pseudoSessions = pseudoSessions_contrast;
        elseif strncmpi(params.trialTypeField,'choice',5) || strncmpi(params.trialTypeField,'movement',5) 
            pseudoSessions = pseudoSessions_choice;
        elseif strncmpi(params.trialTypeField,'feedback',5)
            pseudoSessions = pseudoSessions_feedback;
        elseif strncmpi(params.trialTypeField,'probability',5) || strncmpi(params.trialTypeField,'prior',5)
            try
                pseudoSessions = pseudoSessions_bias;
            catch
                fprintf('Could not find block pseudoSessions, skipping.\n')
                continue
            end
        end
    end

    %overwrite pseudoSessions with new set


    %compute base valid trials
    nTrials = size(allResps,2);
    baseValid = true(1,nTrials);
    baseValid(max(1,nTrials-19):nTrials) = false;           % last 20
    baseValid(any(isnan(allResps),1)) = false;               % any NaN in PETH
    baseValid(all(allResps==0,1)) = false;                   % all-zero trials
    baseValid(badTrials(:)') = false;                       % bad frameQC trials

    if strcmp(params.cv,'odd')
        baseValid(2:2:end) = false;
    elseif strcmp(params.cv,'even')
        baseValid(1:2:end) = false;
    end

    nSess = 1 + nPseudoSessions;
    stats = nan(nROIs, nSess);
    ok    = false(1, nSess);

    sessionsT = cell(1, nSess);
    sessionsT{1} = trialsT;

    % Build pseudo tables once (serial) so parfor has safe indexing
    for k = 1:nPseudoSessions
        sessionsT{k+1} = struct2table(pseudoSessions(k));
    end

    % comparison indices
    if isscalar(params.nComp)
        iComp = 1:params.nComp;
    else
        iComp = params.nComp;
    end
    
    if ~isnan(iComp)
        if ~isempty(params.trialTypeVals)
            if length(params.trialTypeVals)<2*length(iComp)
                fprintf('Length of iComp must be >= 2 * length of selected trial types - skipping this condition.\n')
                continue
            end
        elseif length(unique(trialsT.(params.trialTypeField)))<2*length(iComp)
            fprintf('Not enough unique trial types to make valid comparison, skipping this condition.\n')
            continue
        end
    end


    %pre-allocate error message
    err = strings(1, nSess);

    % pull params fields into locals for parfor robustness
    trialTypeField  = params.trialTypeField;
    condFields      = params.condFields;
    %trialTypeVals   = params.trialTypeVals;
    trialTypeFilter = params.trialTypeFilter;
    nComp           = params.nComp;
    nTrialsMin      = params.nTrialsMin;
    nTrialsToKeep   = params.nTrialsToKeep;
    stat_to_use     = params.stat_to_use;

    fprintf('Parallel-computing ROI test-statistics over %d (pseudo-)sessions..',nSess)

    parfor s = 1:nSess

        try

            % Pick the session trials table deterministically
            sessT = sessionsT{s};

            % trial conditions for this session
            trialCnd = trialSelection(sessT, trialTypeField, trialTypeVals, trialTypeFilter);

            % comp mask
            nCnd = size(trialCnd,1);
            if nCnd>1
                if isscalar(nComp)
                    compMask = sum(trialCnd(1:nComp,:),1) | sum(trialCnd(end-(nComp-1):end,:),1);
                else
                    compMask = sum(trialCnd(nComp,:),1) | sum(trialCnd(1+nCnd-nComp,:),1);
                end
            else
                compMask = trialCnd;
            end

            validMask = baseValid & compMask;

            % trial count criterion
            if nTrialsMin && sum(validMask) < nTrialsMin
                ok(s) = false;
                continue
            end

            % optionally cap trials (keep first N valid)
            if nTrialsToKeep && sum(validMask) > nTrialsToKeep
                tmp = false(size(validMask));
                tmp(find(validMask, nTrialsToKeep, 'first')) = true;  % NOTE: nTrialsToKeep, not nTrialsMin
                validMask = tmp;
            end

            % build labels for the two groups
            validTrialCnd = trialCnd(:, validMask);
            if ~isnan(iComp)
                iValidTrialCnd1 = logical(sum(validTrialCnd(iComp,:),1));
                iValidTrialCnd2 = logical(sum(validTrialCnd(1+numel(trialTypeVals)-iComp,:),1));
            end

            validResps = allResps(:, validMask);

            % compute stats
            switch stat_to_use
                case 'dprime'
                    test_stat = getdprimes(validResps, [iValidTrialCnd2; iValidTrialCnd1]);
                case 'dprime_cv'
                    test_stat = getdprimes(validResps, [iValidTrialCnd2; iValidTrialCnd1], 5);
                case 'auc'
                    test_stat = getAUCs(validResps, [iValidTrialCnd2; iValidTrialCnd1]);
                case 'auc_cv'
                    test_stat = getAUCs(validResps, [iValidTrialCnd2; iValidTrialCnd1], 5);
                case 'meandiff'
                    test_stat = getmeandiffs(validResps, [iValidTrialCnd2; iValidTrialCnd1]);
                case 'ccMean'
                    sessT2 = sessT(validMask,:);
                    sessT2.contrastDiff = sign(sessT2.contrastDiff);
                    test_stat = getccMeans(validResps, sessT2, trialTypeField, condFields);
                case 'ccMeanDiff'
                    sessT2 = sessT(validMask,:);
                    sessT2.contrastDiff = sign(sessT2.contrastDiff);
                    test_stat = getccMeanDiffs(validResps, sessT2, iValidTrialCnd1, trialTypeField, condFields);
                case 'ccu'
                    sessT2 = sessT(validMask,:);
                    sessT2.contrastDiff = sign(sessT2.contrastDiff);
                    test_stat = getccUs(validResps, sessT2, iValidTrialCnd1, trialTypeField, condFields);
                case 'ccMI'
                    sessT2 = sessT(validMask,:);
                    sessT2.contrastDiff = sign(sessT2.contrastDiff);
                    test_stat = getccMIs(validResps, sessT2, iValidTrialCnd1, trialTypeField, condFields);
                otherwise
                    test_stat = nan(nROIs,1);
            end

            stats(:,s) = test_stat(:);
            ok(s) = true;

        catch ME
            err(s) = string(ME.message);
        end

    end
    fprintf('. Done!\n');

    bad = find(err ~= "");
    disp([bad(:), err(bad).'])

    % Real session failed criterion → skip this iParam
    if ~ok(1)
        warning('Test %s produced fewer than %d valid trials in REAL session, skipping this type.\n', stat_names{iParam}, nTrialsMin);
        continue
    end
    % Real session are all nans → skip this iParam
    if all(isnan(stats(:,1)))
        warning('Test %s produced all-NaN in REAL session, skipping this type.', stat_names{iParam});
        continue
    end

    %make sure we take exactly nPseudoSessions
    goodPseudo = find(ok(2:end));
    if numel(goodPseudo) < nPseudoSessions
        warning('Only %d/%d pseudos passed nTrialsMin.', numel(goodPseudo), nPseudoSessions);
    end
    use = goodPseudo(1:min(numel(goodPseudo), nPseudoSessions));
    stats = [stats(:,1), stats(:,1+use)];

    % store the PETH matrix and info
    PETH_struct.entries(iParam).activity_type    = params.activity_type;
    PETH_struct.entries(iParam).evnt             = evnt;
    PETH_struct.entries(iParam).evnt_nm          = evnt_nm;
    %PETH_struct.entries(iParam).twin_all         = params.twin_all;
    %PETH_struct.entries(iParam).twin_bl          = params.twin_bl;
    PETH_struct.entries(iParam).twin_ev          = params.twin_ev;
    %PETH_struct.entries(iParam).T                = T(:)';
    PETH_struct.entries(iParam).badTrials        = badTrials(:);
    PETH_struct.entries(iParam).Resps             = allResps;

    %compute which neurons were significantly responsive to contrast of interest
    stat_pval = nan(1,nROIs);
    for i = 1:nROIs
        stat_pval(i) = 0.01*comp_percentile(stats(i,2:end),stats(i,1));
    end
    if all(isnan(stat_pval))
        warning('Test %s produced all-NaN p-values; skipping column export.', stat_names{iParam});
        continue
    end
    hs = stat_pval>1-params.pthresh/2 | stat_pval<params.pthresh/2; %2-sided

    PETH_struct.out(iParam).stats     = stats;
    PETH_struct.out(iParam).ps        = stat_pval(:);
    PETH_struct.out(iParam).hs        = hs(:);

    statsMat(:, iParam) = stats(:,1);
    psMat(:, iParam) = stat_pval(:);
    hsMat(:, iParam) = hs(:);
    didCompute(iParam) = true;

end



%% save results to each respective FOV folder
if saveToFOV
    fprintf('Saving results..');
    uFov = unique(fov(:))';

    for k = 1:numel(uFov)
        thisFov = uFov(k);

        fovFolder = fullfile(datpath, 'alf', sprintf('FOV_%02d', thisFov));
        if ~exist(fovFolder, 'dir')
            warning('FOV folder not found: %s (skipping)', fovFolder);
            continue
        end

        roiMask = (fov == thisFov);
        if ~any(roiMask)
            warning('No ROIs found for FOV %d (skipping)', thisFov);
            continue
        end

        outPath_p = fullfile(fovFolder, 'mpciROIs.taskTunedP.tsv');
        outPath_stat = fullfile(fovFolder, 'mpciROIs.taskTunedStat.tsv');

        successNames = stat_names(didCompute);

        % if nothing succeeded, nothing to write/append
        if isempty(successNames)
            continue
        end

        ps_fov_success = psMat(roiMask, didCompute);
        stats_fov_success = statsMat(roiMask, didCompute);
        TnewSuccess_p    = array2table(ps_fov_success, 'VariableNames', successNames);
        TnewSuccess_stat = array2table(stats_fov_success, 'VariableNames', successNames);

        if exist(outPath_p, 'file') && exist(outPath_stat, 'file')
            % Read existing
            try
                Told_p = readtable(outPath_p, 'FileType','text', 'Delimiter','\t', 'PreserveVariableNames',true);
                Told_stat = readtable(outPath_stat, 'FileType','text', 'Delimiter','\t', 'PreserveVariableNames',true);
            catch ME
                warning('Failed reading existing %s (%s). Will write new results only.', outPath_p, ME.message);
                Told_p = table();
                Told_stat = table();
            end

            if isempty(Told_p)
                % nothing usable to merge with
                writetable(TnewSuccess_p, outPath_p, 'FileType','text', 'Delimiter','\t');
                writetable(TnewSuccess_stat, outPath_stat, 'FileType','text', 'Delimiter','\t');
                continue
            end

            % row-count safety
            if height(Told_p) ~= height(TnewSuccess_p)
                warning('Row count mismatch in %s (old=%d, new=%d). Will NOT merge; writing new results only.', ...
                    outPath_p, height(Told_p), height(TnewSuccess_p));
                writetable(TnewSuccess_p, outPath_p, 'FileType','text', 'Delimiter','\t');
                writetable(TnewSuccess_stat, outPath_stat, 'FileType','text', 'Delimiter','\t');
                continue
            end

            % column name safety
            if height(Told_p) ~= height(TnewSuccess_p)
                warning('Row count mismatch in %s (old=%d, new=%d). Will NOT merge; writing new results only.', ...
                    outPath_p, height(Told_p), height(TnewSuccess_p));
                writetable(TnewSuccess_p, outPath_p, 'FileType','text', 'Delimiter','\t');
                writetable(TnewSuccess_stat, outPath_stat, 'FileType','text', 'Delimiter','\t');
                continue
            end

            oldNames = Told_p.Properties.VariableNames;
            oldNames_stat = Told_stat.Properties.VariableNames;

            if ~isequaln(oldNames,oldNames_stat)
                warning('Variable name mismatch in %s and %s. Will NOT merge; writing new results only.', ...
                    outPath_p, outPath_stat);
                writetable(TnewSuccess_p, outPath_p, 'FileType','text', 'Delimiter','\t');
                writetable(TnewSuccess_stat, outPath_stat, 'FileType','text', 'Delimiter','\t');
                continue
            end

            newNames = TnewSuccess_p.Properties.VariableNames;

            if overwriteExisting
                % Replace matching columns + append new ones, keep all others
                common = intersect(oldNames, newNames, 'stable');
                add    = setdiff(newNames, oldNames, 'stable');

                Tout_p = Told_p;
                Tout_stat = Told_stat;

                %TEMPORARY HACK to get rid of faulty ccu columns
                %ix_keep = ~startsWith(oldNames, 'ccu');
                %Tout_p = Tout_p(:,ix_keep);
                %Tout_stat = Tout_stat(:,ix_keep);

                % overwrite values for common columns
                for c = 1:numel(common)
                    nm = common{c};
                    Tout_p.(nm) = TnewSuccess_p.(nm);
                    Tout_stat.(nm) = TnewSuccess_stat.(nm);
                end

                % append any genuinely new columns
                if ~isempty(add)
                    Tout_p = [Tout_p, TnewSuccess_p(:, add)];
                    Tout_stat = [Tout_stat, TnewSuccess_stat(:, add)];
                end

            else
                % append only missing columns (your old behavior)
                missingNames = newNames(~ismember(newNames, oldNames));
                if isempty(missingNames)
                    continue
                end

                Tadd_p = TnewSuccess_p(:, missingNames);
                Tadd_stat = TnewSuccess_stat(:, missingNames);
                Tout_p = [Told_p, Tadd_p];
                Tout_stat = [Told_stat, Tadd_stat];

            end

            writetable(Tout_p, outPath_p, 'FileType','text', 'Delimiter','\t');
            writetable(Tout_stat, outPath_stat, 'FileType','text', 'Delimiter','\t');

        else
            % no existing file: write only successful tests
            writetable(TnewSuccess_p, outPath_p, 'FileType','text', 'Delimiter','\t');
            writetable(TnewSuccess_stat, outPath_stat, 'FileType','text', 'Delimiter','\t');
        end
    end

    fprintf('. Done!\n');
end


