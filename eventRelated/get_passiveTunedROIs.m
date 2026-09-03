function PETH_struct = get_passiveTunedROIs(datpath, varargin)
%
% get_passiveTunedROIs
%
% Loads passive stimulus event timings and all ROI traces from all FOVs in
% an IBL mesoscope session.
%
% Tests responses to:
%   valveOn
%   toneOn
%   noiseOn
%
% For each event, activity is averaged in:
%   response window: [ 0.0,  0.4] s
%   baseline window: [-0.4,  0.0] s
%
% Statistical test:
%   A paired permutation test in which response and baseline windows are
%   randomly swapped independently within each trial.
%
% The observed statistic for each ROI is:
%
%   mean(response - baseline)
%
% A positive statistic indicates greater activity during the response
% window than during the baseline window.
%
% Optional inputs:
%
%   'saveToFOV'          false by default
%   'nPermutations'      199 by default
%   'overwriteExisting'  false by default
%   'nTrialsMin'         20 by default
%   'pthresh'            0.05 by default
%   'randomSeed'         [] by default
%
% Example:
%
% PETH_struct = get_passiveTunedROIs( ...
%     'Y:\Subjects\SP080\2026-07-07\002', ...
%     'saveToFOV', true);
%
% Requires:
%   IBL_loadPassiveStimsTable
%   IBL_loadMesoData
%   get_evoked_full
%   readNPY
%
% Samuel Picard / OpenAI, July 2026


%% Options

p = inputParser;

p.addParameter( ...
    'saveToFOV', false, ...
    @(x) islogical(x) || isnumeric(x));

p.addParameter( ...
    'nPermutations', 199, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

p.addParameter( ...
    'overwriteExisting', false, ...
    @(x) islogical(x) || isnumeric(x));

p.addParameter( ...
    'nTrialsMin', 20, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

p.addParameter( ...
    'pthresh', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);

p.addParameter( ...
    'randomSeed', [], ...
    @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x)));

p.parse(varargin{:});

saveToFOV = logical(p.Results.saveToFOV);
nPermutations = round(p.Results.nPermutations);
overwriteExisting = logical(p.Results.overwriteExisting);
nTrialsMin = round(p.Results.nTrialsMin);
pthresh = p.Results.pthresh;
randomSeed = p.Results.randomSeed;

if ~isempty(randomSeed)
    rng(randomSeed);
end


%% Fixed passive-stimulus parameters

eventNames = { ...
    'valveOn', ...
    'toneOn', ...
    'noiseOn'};

statNames = { ...
    'valve', ...
    'tone', ...
    'noise'};

twinResponse = [0, 0.4];
twinBaseline = [-0.4, 0];
twinPETH = [-1, 3];

nEvents = numel(eventNames);


%% Session path information

datpath = char(datpath);
datpath = stripTrailingFilesep(datpath);

[subj, day, sess] = parseSessionPath(datpath);

fprintf('\n%s:\n', datpath);


%% Determine which tests need to run

runMask = true(1, nEvents);

if saveToFOV && ~overwriteExisting

    fovBase = fullfile(datpath, 'alf');
    fovDirs = dir(fullfile(fovBase, 'FOV_*'));

    if ~isempty(fovDirs)
        missingAny = false(1, nEvents);

        for iFov = 1:numel(fovDirs)

            if ~fovDirs(iFov).isdir
                continue
            end

            fovFolder = fullfile( ...
                fovBase, ...
                fovDirs(iFov).name);

            pPath = fullfile( ...
                fovFolder, ...
                'mpciROIs.passiveTunedP.tsv');

            statPath = fullfile( ...
                fovFolder, ...
                'mpciROIs.passiveTunedStat.tsv');

            if exist(pPath, 'file') && exist(statPath, 'file')
                try
                    previousP = readtable( ...
                        pPath, ...
                        'FileType', 'text', ...
                        'Delimiter', '\t', ...
                        'VariableNamingRule', 'preserve');

                    previousStat = readtable( ...
                        statPath, ...
                        'FileType', 'text', ...
                        'Delimiter', '\t', ...
                        'VariableNamingRule', 'preserve');

                    previousNames = intersect( ...
                        previousP.Properties.VariableNames, ...
                        previousStat.Properties.VariableNames);

                    missingHere = ~ismember( ...
                        statNames, ...
                        previousNames);

                catch
                    missingHere = true(1, nEvents);
                end

            else
                missingHere = true(1, nEvents);
            end

            missingAny = missingAny | missingHere;
        end

        runMask = missingAny;
    end

    if ~any(runMask)
        fprintf(['All passive tuning results already exist. ' ...
            'Nothing new to compute.\n']);

        PETH_struct = [];
        return
    end

    fprintf( ...
        'Will compute %d/%d passive stimulus tests.\n', ...
        sum(runMask), ...
        nEvents);
end


%% Load passive stimulus timings

fprintf('Getting passive stimulus timings... ');

try
    passiveT = IBL_loadPassiveStimsTable(datpath);
catch ME
    warning( ...
        'Could not load passive stimulus table:\n%s', ...
        getReport(ME, 'basic'));

    PETH_struct = [];
    return
end

if isempty(passiveT) || height(passiveT) == 0
    fprintf('No passive stimuli found. Skipping this session.\n');

    PETH_struct = [];
    return
end

missingEventFields = eventNames( ...
    ~ismember( ...
    eventNames, ...
    passiveT.Properties.VariableNames));

if ~isempty(missingEventFields)
    warning( ...
        'Passive table is missing event columns: %s', ...
        strjoin(missingEventFields, ', '));
end

fprintf('Done! Found %d passive stimulus rows.\n', height(passiveT));


%% Load neural data

frameQCPath = fullfile( ...
    datpath, ...
    'alf', ...
    'FOV_00', ...
    'mpci.mpciFrameQC.npy');

if exist(frameQCPath, 'file')
    frameQC = readNPY(frameQCPath);

    if any(frameQC ~= 0)
        warning([ ...
            'Bad frames / non-zero frameQC found. ' ...
            'Trials overlapping these frames will be excluded.']);
    end
end

fprintf('Loading 2PI traces... ');

Fall = IBL_loadMesoData( ...
    subj, ...
    day, ...
    sess, ...
    'trace', ...
    'spks', ...
    'fast', ...
    true);

fprintf('Done!\n');

if isempty(Fall.tr) || isempty(Fall.time)
    warning('Incomplete imaging dataset. Skipping this session.');

    PETH_struct = [];
    return
end

if size(Fall.tr, 2) ~= numel(Fall.time)
    warning( ...
        ['Unequal number of frames in ROI activity (%d) ' ...
        'and imaging times (%d). Skipping this session.'], ...
        size(Fall.tr, 2), ...
        numel(Fall.time));

    PETH_struct = [];
    return
end

sig = Fall.tr;
frameTimes = Fall.time(:);

fov = Fall.fov;
idx = Fall.idx;

if isfield(Fall, 'timeshift')
    timeShifts = Fall.timeshift;
else
    timeShifts = [];
end

if isfield(Fall, 'frameQC')
    badframes = find(Fall.frameQC ~= 0);
else
    badframes = [];
end

nROIs = size(sig, 1);


%% Sampling-rate guardrail

minimumWindowWidth = min( ...
    abs(diff(twinResponse)), ...
    abs(diff(twinBaseline)));

[sig, frameTimes, badframes, timeShifts] = ...
    resampleImagingIfNeeded( ...
    sig, ...
    frameTimes, ...
    badframes, ...
    timeShifts, ...
    minimumWindowWidth);

Fs = 1 / median(diff(frameTimes));

%% Passive-period normalization values

onFields = {'valveOn', 'toneOn', 'noiseOn'};
offFields = {'valveOff', 'toneOff', 'noiseOff'};

allOnTimes = [];
allOffTimes = [];

for iField = 1:numel(onFields)

    if ismember(onFields{iField}, passiveT.Properties.VariableNames)
        allOnTimes = [ ...
            allOnTimes; ...
            passiveT.(onFields{iField})(:)]; %#ok<AGROW>
    end

    if ismember(offFields{iField}, passiveT.Properties.VariableNames)
        allOffTimes = [ ...
            allOffTimes; ...
            passiveT.(offFields{iField})(:)]; %#ok<AGROW>
    end
end

allOnTimes = allOnTimes(isfinite(allOnTimes));
allOffTimes = allOffTimes(isfinite(allOffTimes));

if isempty(allOnTimes) || isempty(allOffTimes)
    warning( ...
        ['Could not determine the passive-stimulation period. ' ...
        'Skipping this session.']);

    PETH_struct = [];
    return
end

passiveTimes = [ ...
    min(allOnTimes), ...
    max(allOffTimes)];

passiveFrames = ...
    frameTimes >= passiveTimes(1) & ...
    frameTimes <= passiveTimes(2);

if ~any(passiveFrames)
    warning( ...
        ['No imaging frames were found within the passive-stimulation ' ...
        'period [%.3f, %.3f] s. Skipping this session.'], ...
        passiveTimes(1), passiveTimes(2));

    PETH_struct = [];
    return
end

% Per-ROI normalization limits calculated using activity during the
% passive-stimulation period.
normalizationPrctiles = prctile( ...
    sig(:, passiveFrames), ...
    [20, 99], ...
    2);

normLow = normalizationPrctiles(:, 1);
normHigh = normalizationPrctiles(:, 2);

%% Initialize outputs

PETH_struct = struct();

PETH_struct.meta = struct( ...
    'datpath', datpath, ...
    'subject', subj, ...
    'date', day, ...
    'session', sess, ...
    'Fs', Fs, ...
    'nROIs', nROIs, ...
    'nFrames', size(sig, 2), ...
    'nPermutations', nPermutations, ...
    'nTrialsMin', nTrialsMin, ...
    'pthresh', pthresh, ...
    'responseWindow', twinResponse, ...
    'baselineWindow', twinBaseline);

PETH_struct.roi = struct( ...
    'fov', fov(:), ...
    'idx', idx(:));

PETH_struct.passiveT = passiveT;

PETH_struct.entries = repmat( ...
    struct(), ...
    1, ...
    nEvents);

PETH_struct.out = repmat( ...
    struct(), ...
    1, ...
    nEvents);

statsMat = nan(nROIs, nEvents);
psMat = nan(nROIs, nEvents);
hsMat = false(nROIs, nEvents);
didCompute = false(1, nEvents);


%% Compute responses, full PETHs, and paired permutation tests

% One cell per stimulus.
%
% pethAvgNorm{iEvent}:
%   nROIs x nTimePoints mean normalized PETH
%
% pethTimeValues{iEvent}:
%   1 x nTimePoints vector relative to event onset
%
% pethBadTrials{iEvent}:
%   nTrials x 1 logical vector
pethAvgNorm = cell(1, nEvents);
pethAvgNormOdd = cell(1, nEvents);
pethAvgNormEven = cell(1, nEvents);

pethTimeValues = cell(1, nEvents);
pethBadTrials = cell(1, nEvents);

for iEvent = 1:nEvents

    if ~runMask(iEvent)
        continue
    end

    eventName = eventNames{iEvent};
    statName = statNames{iEvent};

    if ~ismember( ...
            eventName, ...
            passiveT.Properties.VariableNames)

        warning( ...
            'Event column "%s" is not present. Skipping.', ...
            eventName);

        continue
    end

    fprintf( ...
        ['Computing passive responses to %s: ' ...
        'response [%.2f, %.2f] s, ' ...
        'baseline [%.2f, %.2f] s, ' ...
        'full PETH [%.2f, %.2f] s... '], ...
        eventName, ...
        twinResponse(1), ...
        twinResponse(2), ...
        twinBaseline(1), ...
        twinBaseline(2), ...
        twinPETH(1), ...
        twinPETH(2));

    try
        % evokedResponse:
        %   nROIs x nTrials
        %   mean(response window) - mean(baseline window)
        %
        % fullPETH:
        %   nROIs x nTrials x nTimePoints
        %   baseline-subtracted trace for every trial
        %
        % timeValues:
        %   time relative to the locking event
        [evokedResponse, fullPETH, timeValues, badTrials] = ...
            get_evoked_full( ...
            sig, ...
            frameTimes, ...
            badframes, ...
            passiveT, ...
            eventName, ...
            twinResponse, ...
            twinBaseline, ...
            twinPETH, ...
            timeShifts);

    catch ME
        warning( ...
            'Error computing responses for %s:\n%s', ...
            eventName, ...
            getReport(ME, 'basic'));

        continue
    end

    fprintf('Done!\n');

    badTrials = badTrials(:)';

    % Exclude globally invalid trials. Additional NaN handling is done
    % independently for each ROI inside pairedWindowSwapTest.
    globallyValidTrials = ~badTrials;

    if sum(globallyValidTrials) < nTrialsMin
        warning( ...
            ['Only %d globally valid %s trials found; ' ...
            'at least %d are required. Skipping.'], ...
            sum(globallyValidTrials), ...
            eventName, ...
            nTrialsMin);

        continue
    end


    %% Compute all-trial, odd-trial, and even-trial mean PETHs

    % Ensure globally invalid trials do not contribute.
    fullPETH(:, badTrials, :) = NaN;

    nTrials = size(fullPETH, 2);

    oddTrialMask = false(1, nTrials);
    oddTrialMask(1:2:end) = true;

    evenTrialMask = false(1, nTrials);
    evenTrialMask(2:2:end) = true;

    % Also exclude trials marked bad.
    oddTrialMask = oddTrialMask & ~badTrials;
    evenTrialMask = evenTrialMask & ~badTrials;


    % All valid trials
    meanPETH = meanPETHAcrossTrials( ...
        fullPETH, ...
        ~badTrials);


    % Odd-numbered valid trials
    meanPETHOdd = meanPETHAcrossTrials( ...
        fullPETH, ...
        oddTrialMask);


    % Even-numbered valid trials
    meanPETHEven = meanPETHAcrossTrials( ...
        fullPETH, ...
        evenTrialMask);


    % Apply the same passive-period normalization bounds to all versions.
    meanPETHNorm = normalize_minmax( ...
        meanPETH, ...
        normLow, ...
        normHigh);

    meanPETHNormOdd = normalize_minmax( ...
        meanPETHOdd, ...
        normLow, ...
        normHigh);

    meanPETHNormEven = normalize_minmax( ...
        meanPETHEven, ...
        normLow, ...
        normHigh);


    pethAvgNorm{iEvent} = single(meanPETHNorm);
    pethAvgNormOdd{iEvent} = single(meanPETHNormOdd);
    pethAvgNormEven{iEvent} = single(meanPETHNormEven);

    pethTimeValues{iEvent} = single(timeValues(:)');
    pethBadTrials{iEvent} = badTrials(:);

    %% Paired response-versus-baseline permutation test

    % evokedResponse is already:
    %
    %   response window - baseline window
    %
    % Swapping response and baseline within a trial is therefore
    % equivalent to multiplying the paired difference by -1.
    fprintf( ...
        'Running %d paired window-swap permutations for %s... ', ...
        nPermutations, ...
        eventName);

    %normalize first
    evokedResponseNorm = normalize_minmax(evokedResponse, normLow, normHigh);

    [observedStat, nullStats, pValues, significant, nValidTrials] = ...
        pairedWindowSwapTest( ...
        evokedResponseNorm(:, globallyValidTrials), ...
        nPermutations, ...
        nTrialsMin, ...
        pthresh);

    fprintf('Done!\n');

    if all(isnan(observedStat))
        warning( ...
            'All observed statistics were NaN for %s. Skipping.', ...
            eventName);

        continue
    end


    %% Store event-level data and provenance

    PETH_struct.entries(iEvent).name = statName;
    PETH_struct.entries(iEvent).evnt = eventName;

    PETH_struct.entries(iEvent).twin_ev = twinResponse;
    PETH_struct.entries(iEvent).twin_bl = twinBaseline;
    PETH_struct.entries(iEvent).twin_all = twinPETH;

    PETH_struct.entries(iEvent).badTrials = badTrials(:);

    % Baseline-subtracted response-window value for every ROI and trial.
    PETH_struct.entries(iEvent).evokedResponse = evokedResponse;
    PETH_struct.entries(iEvent).evokedResponseNorm = evokedResponseNorm;

    % Full baseline-subtracted ROI x trial x time PETH.
    PETH_struct.entries(iEvent).fullPETH = ...
        fullPETH;

    % Mean normalized ROI x time PETH.
    PETH_struct.entries(iEvent).PETHavgNorm = ...
        pethAvgNorm{iEvent};

    PETH_struct.entries(iEvent).PETHavgNorm_odd = ...
        pethAvgNormOdd{iEvent};

    PETH_struct.entries(iEvent).PETHavgNorm_even = ...
        pethAvgNormEven{iEvent};

    PETH_struct.entries(iEvent).timeValues = ...
        pethTimeValues{iEvent};

    PETH_struct.out(iEvent).stats = observedStat(:);
    PETH_struct.out(iEvent).nullStats = nullStats;
    PETH_struct.out(iEvent).ps = pValues(:);
    PETH_struct.out(iEvent).hs = significant(:);
    PETH_struct.out(iEvent).nValidTrials = nValidTrials(:);

    statsMat(:, iEvent) = observedStat(:);
    psMat(:, iEvent) = pValues(:);
    hsMat(:, iEvent) = significant(:);

    didCompute(iEvent) = true;
end


%% Save results in each FOV folder

if saveToFOV

    fprintf('Saving passive tuning results... ');

    successfulNames = statNames(didCompute);

    if isempty(successfulNames)
        fprintf('No tests completed successfully. Nothing to save.\n');
        return
    end

    uniqueFOVs = unique(fov(:))';

    for iFov = 1:numel(uniqueFOVs)

        thisFOV = uniqueFOVs(iFov);

        fovFolder = fullfile( ...
            datpath, ...
            'alf', ...
            sprintf('FOV_%02d', thisFOV));

        if ~exist(fovFolder, 'dir')
            warning( ...
                'FOV folder not found: %s. Skipping.', ...
                fovFolder);

            continue
        end

        roiMask = fov == thisFOV;

        if ~any(roiMask)
            continue
        end

        %% Save all-trial, odd-trial, and even-trial passive PETHs

        for iEvent = 1:nEvents

            if ~didCompute(iEvent)
                continue
            end

            if isempty(pethAvgNorm{iEvent}) || ...
                    isempty(pethAvgNormOdd{iEvent}) || ...
                    isempty(pethAvgNormEven{iEvent}) || ...
                    isempty(pethTimeValues{iEvent})
                continue
            end

            eventName = eventNames{iEvent};

            pethAllForFOV = pethAvgNorm{iEvent}(roiMask, :);
            pethOddForFOV = pethAvgNormOdd{iEvent}(roiMask, :);
            pethEvenForFOV = pethAvgNormEven{iEvent}(roiMask, :);

            timeValues = pethTimeValues{iEvent};


            allPath = fullfile( ...
                fovFolder, ...
                sprintf( ...
                'mpciROIs.PETHavgNorm_%s_passive.npy', ...
                eventName));

            oddPath = fullfile( ...
                fovFolder, ...
                sprintf( ...
                'mpciROIs.PETHavgNorm_%s_passive_odd.npy', ...
                eventName));

            evenPath = fullfile( ...
                fovFolder, ...
                sprintf( ...
                'mpciROIs.PETHavgNorm_%s_passive_even.npy', ...
                eventName));

            timePath = fullfile( ...
                fovFolder, ...
                sprintf( ...
                'PETHavgNorm_%s_passive.timeValues.npy', ...
                eventName));


            writeNPY(single(pethAllForFOV), allPath);
            writeNPY(single(pethOddForFOV), oddPath);
            writeNPY(single(pethEvenForFOV), evenPath);
            writeNPY(single(timeValues(:)), timePath);

        end

        pPath = fullfile( ...
            fovFolder, ...
            'mpciROIs.passiveTunedP.tsv');

        statPath = fullfile( ...
            fovFolder, ...
            'mpciROIs.passiveTunedStat.tsv');

        % hPath = fullfile( ...
        %     fovFolder, ...
        %     'mpciROIs.passiveTunedH.tsv');

        newP = array2table( ...
            psMat(roiMask, didCompute), ...
            'VariableNames', successfulNames);

        newStat = array2table( ...
            statsMat(roiMask, didCompute), ...
            'VariableNames', successfulNames);

        % newH = array2table( ...
        %     hsMat(roiMask, didCompute), ...
        %     'VariableNames', successfulNames);

        mergeOrWriteResultTable( ...
            pPath, ...
            newP, ...
            overwriteExisting);

        mergeOrWriteResultTable( ...
            statPath, ...
            newStat, ...
            overwriteExisting);

        % mergeOrWriteResultTable( ...
        %     hPath, ...
        %     newH, ...
        %     overwriteExisting);
    end

    fprintf('Done!\n');
end

end



function [sig, frameTimes, badframes, timeShifts] = ...
    resampleImagingIfNeeded( ...
    sig, ...
    frameTimes, ...
    badframes, ...
    timeShifts, ...
    minimumWindowWidth)

frameTimes = frameTimes(:);

dt = median(diff(frameTimes));

if dt <= minimumWindowWidth
    return
end

% Follow the same guardrail approach as get_taskTunedROIs: choose a new
% timestep slightly below the original sampling interval, rounded down to
% a multiple of 0.05 seconds.
newDt = floor(20 * (0.99 * dt)) / 20;

% Ensure that resampling actually increases temporal resolution.
if newDt <= 0 || newDt >= dt
    newDt = 0.99 * dt;
end

warning( ...
    ['Sampling interval %.4f s is larger than or equal to the ' ...
    'shortest analysis window %.4f s. Resampling to %.4f s ' ...
    'without interpolating across acquisition gaps.'], ...
    dt, ...
    minimumWindowWidth, ...
    newDt);

originalFrameTimes = frameTimes;
originalBadframes = badframes;

gapThreshold = 3 * dt;
gapIndices = find(diff(originalFrameTimes) > gapThreshold);

segmentStarts = [1; gapIndices(:) + 1];
segmentEnds = [gapIndices(:); numel(originalFrameTimes)];

newFrameTimeCells = cell(numel(segmentStarts), 1);
newSignalCells = cell(numel(segmentStarts), 1);
newTimeShiftCells = cell(numel(segmentStarts), 1);

hasTimeShifts = ~isempty(timeShifts);

for iSegment = 1:numel(segmentStarts)

    segmentIndices = ...
        segmentStarts(iSegment):segmentEnds(iSegment);

    if numel(segmentIndices) < 2
        continue
    end

    segmentTimes = originalFrameTimes(segmentIndices);
    segmentSignal = sig(:, segmentIndices);

    newSegmentTimes = ...
        (segmentTimes(1):newDt:segmentTimes(end))';

    if isempty(newSegmentTimes)
        continue
    end

    newFrameTimeCells{iSegment} = newSegmentTimes;

    newSignalCells{iSegment} = interp1( ...
        segmentTimes, ...
        segmentSignal', ...
        newSegmentTimes, ...
        'linear')';

    if hasTimeShifts

        if isvector(timeShifts)
            originalShift = timeShifts(:);

            newTimeShiftCells{iSegment} = interp1( ...
                segmentTimes, ...
                originalShift(segmentIndices), ...
                newSegmentTimes, ...
                'linear');

        else
            segmentShift = timeShifts(:, segmentIndices);

            newTimeShiftCells{iSegment} = interp1( ...
                segmentTimes, ...
                segmentShift', ...
                newSegmentTimes, ...
                'linear')';
        end
    end
end

keepSegments = ~cellfun( ...
    @isempty, ...
    newFrameTimeCells);

frameTimes = vertcat( ...
    newFrameTimeCells{keepSegments});

sig = cat( ...
    2, ...
    newSignalCells{keepSegments});

if hasTimeShifts
    if isvector(timeShifts)
        timeShifts = vertcat( ...
            newTimeShiftCells{keepSegments});
    else
        timeShifts = cat( ...
            2, ...
            newTimeShiftCells{keepSegments});
    end
end

badframes = [];

% Remap previously bad imaging frames.
if ~isempty(originalBadframes)

    originalBadTimes = ...
        originalFrameTimes(originalBadframes);

    nearestIndices = interp1( ...
        frameTimes, ...
        1:numel(frameTimes), ...
        originalBadTimes, ...
        'nearest', ...
        nan);

    nearestIndices = nearestIndices(isfinite(nearestIndices));
    badframes = unique(round(nearestIndices(:)))';
end

% Mark frames around original acquisition gaps as bad.
gapBuffer = newDt;

for iGap = 1:numel(gapIndices)

    gapStart = originalFrameTimes(gapIndices(iGap));
    gapEnd = originalFrameTimes(gapIndices(iGap) + 1);

    nearGap = ...
        frameTimes >= gapStart - gapBuffer & ...
        frameTimes <= gapEnd + gapBuffer;

    badframes = unique( ...
        [badframes, find(nearGap)']);
end

fprintf( ...
    ['Resampled from %d to %d frames across ' ...
    '%d continuous imaging segments.\n'], ...
    numel(originalFrameTimes), ...
    numel(frameTimes), ...
    sum(keepSegments));

end


function mergeOrWriteResultTable( ...
    outputPath, ...
    newTable, ...
    overwriteExisting)

if ~exist(outputPath, 'file')

    writetable( ...
        newTable, ...
        outputPath, ...
        'FileType', 'text', ...
        'Delimiter', '\t');

    return
end

try
    oldTable = readtable( ...
        outputPath, ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'VariableNamingRule', 'preserve');

catch ME
    warning( ...
        ['Could not read existing result file %s: %s. ' ...
        'Writing the new results only.'], ...
        outputPath, ...
        ME.message);

    writetable( ...
        newTable, ...
        outputPath, ...
        'FileType', 'text', ...
        'Delimiter', '\t');

    return
end

if height(oldTable) ~= height(newTable)
    warning( ...
        ['Row-count mismatch in %s: old = %d, new = %d. ' ...
        'Writing the new results only.'], ...
        outputPath, ...
        height(oldTable), ...
        height(newTable));

    writetable( ...
        newTable, ...
        outputPath, ...
        'FileType', 'text', ...
        'Delimiter', '\t');

    return
end

oldNames = oldTable.Properties.VariableNames;
newNames = newTable.Properties.VariableNames;

if overwriteExisting

    outputTable = oldTable;

    commonNames = intersect( ...
        oldNames, ...
        newNames, ...
        'stable');

    additionalNames = setdiff( ...
        newNames, ...
        oldNames, ...
        'stable');

    for iName = 1:numel(commonNames)
        name = commonNames{iName};
        outputTable.(name) = newTable.(name);
    end

    if ~isempty(additionalNames)
        outputTable = [ ...
            outputTable, ...
            newTable(:, additionalNames)];
    end

else

    missingNames = newNames( ...
        ~ismember(newNames, oldNames));

    if isempty(missingNames)
        return
    end

    outputTable = [ ...
        oldTable, ...
        newTable(:, missingNames)];
end

writetable( ...
    outputTable, ...
    outputPath, ...
    'FileType', 'text', ...
    'Delimiter', '\t');

end


function datpath = stripTrailingFilesep(datpath)

while ~isempty(datpath) && ...
        (datpath(end) == '/' || datpath(end) == '\')

    datpath(end) = [];
end

end


function [subject, dateString, session] = ...
    parseSessionPath(datpath)

pathParts = split( ...
    string(datpath), ...
    ["/", "\"]);

pathParts(pathParts == "") = [];

if numel(pathParts) < 3
    error( ...
        'get_passiveTunedROIs:InvalidSessionPath', ...
        ['datpath must end with the hierarchy ' ...
        '<subject>\<date>\<session>.']);
end

subject = char(pathParts(end - 2));
dateString = char(pathParts(end - 1));
session = char(pathParts(end));

end

function meanPETH = meanPETHAcrossTrials(fullPETH, trialMask)
%MEANPETHACROSSTRIALS Average an ROI x trial x time PETH over selected trials.
%
% Returns an ROI x time matrix.

nROIs = size(fullPETH, 1);
nTime = size(fullPETH, 3);

trialMask = logical(trialMask(:)');

if numel(trialMask) ~= size(fullPETH, 2)
    error( ...
        'meanPETHAcrossTrials:MaskSizeMismatch', ...
        'trialMask length must equal the number of trials.');
end

if ~any(trialMask)
    meanPETH = nan(nROIs, nTime);
    return
end

meanPETH = mean( ...
    fullPETH(:, trialMask, :), ...
    2, ...
    'omitnan');

meanPETH = reshape( ...
    meanPETH, ...
    nROIs, ...
    nTime);

end