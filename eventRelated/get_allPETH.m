function PETH_struct = get_allPETH(datpath, varargin)
%get_allPETH computes a set of PETHs for all FOVs from a given session
%
% each PETH is a matrix of nROIs x nConds x nTimepoints
%
% Samuel Picard

%% options
p = inputParser;
p.addParameter('saveToFOV', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('cv', false, @(x)islogical(x) || isnumeric(x));
%p.addParameter('saveFilename', 'PETH_struct.mat', @(s)ischar(s) || isstring(s));
p.parse(varargin{:});
saveToFOV = logical(p.Results.saveToFOV);
cv = logical(p.Results.cv);
%saveFilename = char(p.Results.saveFilename);

%define the event types, windows etc to use for the PETHs
params_all = struct( ...
    'activity_type',   'deconv', ...
    'trialTypeField',  {'contrastDiff','choice','feedbackType','probabilityLeft'}, ...
    'trialTypeVals',   {[-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1] [-1,1], [-1,1], [0.2,0.8]}, ...
    ...%'trialTypeVals',   {[-1,0,1] [-1,1], [-1,1], [0.2,0.8]}, ...
    'trialTypeFilter', {'contrastDiff~=0', '', 'contrastDiff~=0', ''}, ...
    'cndFields',       {{'choice'},{'contrastDiff'},{'choice'},{'contrastDiff'}},...
    ...%'cndFields',       {'','','',''}, ...
    'evnt',            {'stimOn','choiceMovement','feedback','stimOn'}, ...
    'twin_all',        {[-1,3],[-2.5,1.5],[-1.5,2.5],[-2,2]}, ...
    'twin_bl',         'none', ...
    'twin_ev',         {[0,0.4],[-0.2,0.2],[0,0.4],[-0.5,-0.1]}, ...
    'nTrialsMin',      50, ...
    'nTrialsToKeep',   false, ...
    'minTrialsPerCond',10, ...
    'minTrialsPerCombo',2 ...
    );

% params_all = struct( ...
%     'activity_type',   'deconv', ...
%     'trialTypeField',  {'contrastDiff','choice','feedbackType'}, ...
%     ...%'trialTypeVals',   {[-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1] [-1,1], [-1,1]}, ...
%     'trialTypeVals',   {[-1,1] [-1,1], [-1,1]}, ...
%     'trialTypeFilter', {'contrastDiff~=0', '', 'contrastDiff~=0'}, ...
%     'cndFields',       {{'choice'},{'contrastDiff'},{'choice'}},...
%     ...'cndFields',       {'','','',''}, ...
%     'evnt',            {'stimOn','choiceMovement','feedback'}, ...
%     'twin_all',        {[-1,3],[-2.5,1.5],[-1.5,2.5]}, ...
%     'twin_bl',         'none', ...
%     'twin_ev',         {[0,0.4],[-0.2,0.2],[0,0.4]}, ...
%     'nTrialsMin',      50, ...
%     'nTrialsToKeep',   false, ...
%     'minTrialsPerCond',10, ...
%     'minTrialsPerCombo',2 ...
%     );

% params_all = struct( ...
%     'activity_type',   'deconv', ...
%     'trialTypeField',  {'choice'}, ...
%     'trialTypeVals',   {[-1,1]}, ...
%     'trialTypeFilter', {'contrastDiff~=0'}, ...
%     'cndFields',       {''}, ...
%     'evnt',            {'firstMovement'}, ...
%     'twin_all',        {[-1,3]}, ...
%     'twin_bl',         {[-0.5,-0.1]}, ...
%     'twin_ev',         {[0,0.4]}, ...
%     'nTrialsMin',      50, ...
%     'nTrialsToKeep',   false, ...
%     'minTrialsPerCond',10, ...
%     'minTrialsPerCombo',2 ...
%     );

% params_all = struct( ...
%     'activity_type',   'deconv', ...
%     'trialTypeField',  {'feedbackType'}, ...
%     'trialTypeVals',   {[-1,1]}, ...
%     'trialTypeFilter', {'contrastDiff~=0'}, ...
%     'cndFields',       {'choice'}, ...
%     'evnt',            {'feedback'}, ...
%     'twin_all',        {[-1.5,2.5]}, ...
%     'twin_bl',         'none', ...
%     'twin_ev',         {[0,0.4]}, ...
%     'nTrialsMin',      50, ...
%     'nTrialsToKeep',   false, ...
%     'minTrialsPerCond',10, ...
%     'minTrialsPerCombo',2 ...
%     );

%path logic
splitPath = split(datpath, filesep);
subj = splitPath{end-2};
dt = splitPath{end-1};
sess = splitPath{end};

fprintf('\n%s:\n ',datpath);

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
elseif size(trialsT,1)<params_all(1).nTrialsMin
    fprintf('fewer than %d total trials, skipping this session!\n', params_all(1).nTrialsMin);
    PETH_struct = [];
    return
end
fprintf('Done!\n');

%% load neural data

%load frame QC, check if there are any non-zero, and throw warning if so
frameQC = readNPY(fullfile(datpath,'alf','FOV_00','mpci.mpciFrameQC.npy'));
if any(frameQC~=0)
    warning('Badframes / non-zero frameQC found. We will ignore trials in those periods.\n');
end

%load data from all FOVs into a single struct
fprintf('Loading 2PI traces...');
Fall = IBL_loadMesoData(subj,dt,sess,'trace','spks','fast',true);

if isempty(Fall.tr) || isempty(Fall.time)
    fprintf('Incomplete datasets, skipping this session!\n');
    PETH_struct = [];
    return
end

sig = Fall.tr;            % deconvolved activity
frameTimes = Fall.time;     % relies on existing mpci.times.npy
Fs = 1/median(diff(Fall.time)); % calculate frameRate from time vector
fov = Fall.fov;             % per-ROI fov label (typically 0..N-1)
idx = Fall.idx;             % per-ROI local index within FOV (if provided)

%additionally derive the badframes from frameQC
badframes = find(Fall.frameQC~=0); % for now, assume all nonzero frameQC is bad

%fprintf(' Done!\n');

%% compute mean PETHs

%initialize PETH struct (session-level metadata + one entry per param)
PETH_struct = struct();
PETH_struct.meta = struct( ...
    'datpath', datpath, ...
    'subject', subj, ...
    'date',    dt, ...
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

for iParam = 1:length(params_all)

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
    if isempty(params.trialTypeVals)
        vals = trialsT.(params.trialTypeField);
        trialTypeVals = unique(vals(~isnan(vals)));
    else
        trialTypeVals = params.trialTypeVals;
    end

    % compute PETHs (raw and time-averaged event triggered responses)
    fprintf('Computing PETHs relative to %s... ', evnt);
    try
        [~, allPETH_raw, T, badTrials] = get_evoked_full( ...
            sig, frameTimes, badframes, trialsT, evnt_nm, ...
            params.twin_ev, params.twin_bl, params.twin_all);
    catch ME
        warning('Error computing PETHs, skipping this type!\n%s\n', getReport(ME, 'basic'));
        %PETH_struct = [];
        continue;
    end

    fprintf('Normalizing and computing averages... ')

    % normalize PETHs (per-ROI) using the task period only
    taskTimes = [trialsT(1,:).intervals_0, trialsT(end,:).intervals_1];
    taskFrames = frameTimes>taskTimes(1) & frameTimes<taskTimes(end);
    prctiles = prctile(sig(:,taskFrames), [20, 99], 2);
    allPETH_raw_norm = normalize_minmax(allPETH_raw, prctiles(:,1), prctiles(:,2));

    %take out badTrials and apply trialTypeFilter
    trialKeep = ~badTrials(:);
    filterMask = evalTrialTypeFilter(trialsT, params.trialTypeFilter);
    trialKeep = trialKeep & filterMask(:);

    allPETH_raw_norm_qc = allPETH_raw_norm(:,trialKeep,:);
    trialsT_qc = trialsT(trialKeep,:);

    if isempty(trialsT_qc)
        warning('No trials left after badTrials and trialTypeFilter for %s. Skipping this type.', evnt);
        continue
    end

    %% Recode trial conditions where required

    effectiveTrialTypeVals = trialTypeVals;

    % Collapse all non-zero contrasts into stimulus side:
    %   negative contrast -> -1
    %   positive contrast -> +1
    %
    % trialTypeFilter='contrastDiff~=0' has already removed zero contrasts.
    if strcmp(params.trialTypeField, 'contrastDiff')

        trialsT_qc.contrastDiff = sign( ...
            trialsT_qc.contrastDiff);

        effectiveTrialTypeVals = [-1, 1];
    end


    %% Count trials for the effective conditions

    nConds = numel(effectiveTrialTypeVals);
    nTrialsPerCond = zeros(1, nConds);

    try
        ttvals = trialsT_qc.(params.trialTypeField);

        for c = 1:nConds
            nTrialsPerCond(c) = sum( ...
                ttvals == effectiveTrialTypeVals(c));
        end

    catch ME
        warning( ...
            'Could not count trials per condition for %s: %s', ...
            params.trialTypeField, ...
            ME.message);

        nTrialsPerCond(:) = NaN;
    end

    insufficientConditions = ...
        nTrialsPerCond < params.minTrialsPerCond;


    %% Compute all-trial mean PETH

    PETH = getMeanResps_allROIs( ...
        allPETH_raw_norm_qc, ...
        trialsT_qc, ...
        params.trialTypeField, ...
        effectiveTrialTypeVals, ...
        params.cndFields, ...
        params.minTrialsPerCombo);

    PETH(:, insufficientConditions, :) = NaN;


    %% Compute odd/even mean PETHs

    if cv

        oddIdx = 1:2:size(trialsT_qc, 1);
        evenIdx = 2:2:size(trialsT_qc, 1);

        PETH_odd = getMeanResps_allROIs( ...
            allPETH_raw_norm_qc(:, oddIdx, :), ...
            trialsT_qc(oddIdx, :), ...
            params.trialTypeField, ...
            effectiveTrialTypeVals, ...
            params.cndFields, ...
            params.minTrialsPerCombo);

        PETH_even = getMeanResps_allROIs( ...
            allPETH_raw_norm_qc(:, evenIdx, :), ...
            trialsT_qc(evenIdx, :), ...
            params.trialTypeField, ...
            effectiveTrialTypeVals, ...
            params.cndFields, ...
            params.minTrialsPerCombo);

        % Use the full retained-trial count to decide whether a stimulus side
        % exists in sufficient numbers. Thus a side is removed only when its
        % total number of trials is below minTrialsPerCond.
        PETH_odd(:, insufficientConditions, :) = NaN;
        PETH_even(:, insufficientConditions, :) = NaN;
    end

    % store params in the PETH_struct
    PETH_struct.entries(iParam).activity_type    = params.activity_type;
    PETH_struct.entries(iParam).trialTypeField   = params.trialTypeField;
    PETH_struct.entries(iParam).trialTypeVals    = effectiveTrialTypeVals;
    PETH_struct.entries(iParam).originalTrialTypeVals = trialTypeVals;
    PETH_struct.entries(iParam).trialTypeFilter  = params.trialTypeFilter;
    PETH_struct.entries(iParam).evnt             = evnt;
    PETH_struct.entries(iParam).evnt_nm          = evnt_nm;
    PETH_struct.entries(iParam).twin_all         = params.twin_all;
    PETH_struct.entries(iParam).twin_bl          = params.twin_bl;
    PETH_struct.entries(iParam).twin_ev          = params.twin_ev;
    PETH_struct.entries(iParam).T                = T(:)';
    PETH_struct.entries(iParam).badTrials        = badTrials(:);
    PETH_struct.entries(iParam).nTrialsPerCond   = nTrialsPerCond;
    PETH_struct.entries(iParam).minTrialsPerCond = params.minTrialsPerCond;

    % store the mean PETH (nROIs x nConds x nTimepoints)
    PETH_struct.entries(iParam).PETH             = PETH;
    if cv
        PETH_struct.entries(iParam).PETH_odd         = PETH_odd;
        PETH_struct.entries(iParam).PETH_even        = PETH_even;
    end

    fprintf('Done!\n');
end

%% save these results to each respective FOV folder
if saveToFOV
    uFov = unique(fov(:))';

    % helper to make safe filename tokens
    makeSafe = @(s) regexprep(char(s), '[^A-Za-z0-9_\-]', '');

    for k = 1:numel(uFov)
        thisFov = uFov(k);

        % typical IBL: alf/FOV_00, alf/FOV_01, ...
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

        % write one file per PETH type
        for iParam = 1:numel(PETH_struct.entries)
            entry = PETH_struct.entries(iParam);

            if ~isfield(entry, 'PETH') || isempty(entry.PETH)
                continue
            end
            if ~isfield(entry, 'T') || isempty(entry.T)
                warning('Missing T for entry %d (skipping metadata save)', iParam);
            end
            if ~isfield(entry, 'trialTypeVals') || isempty(entry.trialTypeVals)
                warning('Missing trialTypeVals for entry %d (skipping metadata save)', iParam);
            end

            % subset to ROIs in this FOV
            PETH_fov = entry.PETH(roiMask, :, :);  % nROIs x nConds x nTimepoints
            PETH_odd_fov = entry.PETH_odd(roiMask, :, :);
            PETH_even_fov = entry.PETH_even(roiMask, :, :);

            % build filename pattern
            evntTok  = makeSafe(entry.evnt);
            ttTok    = makeSafe(entry.trialTypeField);
            fName = sprintf('PETHavgNorm_%s_%s', evntTok, ttTok);

            % PETH file and metadata files (time + condition values)
            outPath = fullfile(fovFolder, ['mpciROIs.', fName]);
            outPathT = fullfile(fovFolder, [fName '.timeValues.npy']);
            outPathC = fullfile(fovFolder, [fName '.conditionValues.npy']);

            try
                writeNPY(PETH_fov, [outPath, '.npy']);
                if cv
                    writeNPY(PETH_odd_fov, [outPath, '_odd.npy']);
                    writeNPY(PETH_even_fov, [outPath, '_even.npy']);
                end
                if isfield(entry, 'T') && ~isempty(entry.T)
                    writeNPY(entry.T(:), outPathT);                 % nTimepoints
                end
                if isfield(entry, 'trialTypeVals') && ~isempty(entry.trialTypeVals)
                    writeNPY(entry.trialTypeVals(:), outPathC);     % nConditions
                end
            catch ME
                warning('Failed writing %s\n%s', outPath, getReport(ME,'basic'));
            end
        end
    end
end

