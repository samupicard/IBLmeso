function out = get_passiveMovieCorrROIs(datpath, varargin)
%GET_PASSIVEMOVIECORRROIS Compute passive-movie cross-repeat reliability.
%
% out = get_passiveMovieCorrROIs(datpath)
%
% For each FOV in:
%   <datpath>/alf/FOV_XX/
%
% this function loads:
%   mpci.ROIActivityDeconvolved.npy
%
% It also loads:
%   mpci.times.npy
%   _sp_video.times.npy
%
% and computes, for each ROI, the Pearson correlation between its activity
% during passive-movie repeat 1 and repeat 2.
%
% Results are saved in each FOV folder as:
%   mpciROIs.passiveMovieCorr.npy
%
% The output vector has one value per ROI, in the same column order as
% mpci.ROIActivityDeconvolved.npy.
%
% Assumptions
% -----------
% - Traces are arranged nFrames x nROIs. If the first dimension does not
%   match frameTimes but the second does, the matrix is transposed.
% - _sp_video.times.npy is nVideoFrames x nRepeats.
% - At least two repeats are present.
%
% Example
% -------
% out = get_passiveMovieCorrROIs( ...
%     'Y:\Subjects\SP072\2025-09-02\002', ...
%     'overwriteExisting', true);
%
% Samuel Picard / OpenAI, July 2026


%% Options

p = inputParser;

p.addParameter( ...
    'TraceFile', ...
    'mpci.ROIActivityDeconvolved.npy', ...
    @(x) ischar(x) || isstring(x));

p.addParameter( ...
    'FrameTimesFile', ...
    'mpci.times.npy', ...
    @(x) ischar(x) || isstring(x));

p.addParameter( ...
    'VideoTimesFile', ...
    '_sp_video.times.npy', ...
    @(x) ischar(x) || isstring(x));

p.addParameter( ...
    'OutputFile', ...
    'mpciROIs.passiveMovieCorr.npy', ...
    @(x) ischar(x) || isstring(x));

p.addParameter( ...
    'FOVs', ...
    [], ...
    @(x) isempty(x) || isnumeric(x));

p.addParameter( ...
    'DurationToleranceFrames', ...
    1, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);

p.addParameter( ...
    'overwriteExisting', ...
    false, ...
    @(x) islogical(x) || isnumeric(x));

p.addParameter( ...
    'MakePlots', ...
    false, ...
    @(x) islogical(x) || isnumeric(x));

p.parse(varargin{:});
opts = p.Results;

opts.TraceFile = char(opts.TraceFile);
opts.FrameTimesFile = char(opts.FrameTimesFile);
opts.VideoTimesFile = char(opts.VideoTimesFile);
opts.OutputFile = char(opts.OutputFile);
opts.overwriteExisting = logical(opts.overwriteExisting);
opts.MakePlots = logical(opts.MakePlots);


%% Session-level file discovery

fprintf('\n%s:\n', datpath);

alfPath = fullfile(datpath, 'alf');

if ~isfolder(alfPath)
    error('ALF folder not found: %s', alfPath);
end

fovDirs = dir(fullfile(alfPath, 'FOV_*'));
fovDirs = fovDirs([fovDirs.isdir]);

if isempty(fovDirs)
    warning('No FOV folders found in %s.', alfPath);
    out = struct([]);
    return
end

availableFOVs = nan(1, numel(fovDirs));

for iFOV = 1:numel(fovDirs)
    availableFOVs(iFOV) = sscanf(fovDirs(iFOV).name, 'FOV_%d');
end

availableFOVs = sort(availableFOVs(isfinite(availableFOVs)));

if isempty(opts.FOVs)
    fovs = availableFOVs;
else
    fovs = intersect(opts.FOVs(:)', availableFOVs, 'stable');

    missingFOVs = setdiff(opts.FOVs(:)', availableFOVs);

    if ~isempty(missingFOVs)
        warning( ...
            'Requested FOVs not found in %s: %s', ...
            datpath, mat2str(missingFOVs));
    end
end

if isempty(fovs)
    warning('No requested FOVs were found. Nothing to compute.');
    out = struct([]);
    return
end


%% Load neural frame times

frameTimesPath = findFrameTimesFile(alfPath, fovs, opts.FrameTimesFile);

if isempty(frameTimesPath)
    warning('Could not find %s in any requested FOV.', opts.FrameTimesFile);
    out = struct([]);
    return
end

frameTimes = double(readNPY(frameTimesPath));
frameTimes = frameTimes(:);

if isempty(frameTimes) || any(~isfinite(frameTimes))
    error('Invalid neural frame times in %s.', frameTimesPath);
end


%% Load passive-movie video times

videoTimesPath = findFileRecursive(alfPath, opts.VideoTimesFile);

if isempty(videoTimesPath)
    warning('Could not find %s below %s.', opts.VideoTimesFile, alfPath);
    out = struct([]);
    return
end

videoTimes = double(readNPY(videoTimesPath));

if isvector(videoTimes)
    warning( ...
        '%s contains only one vector; at least two repeats are required.', ...
        videoTimesPath);
    out = struct([]);
    return
end

if size(videoTimes, 2) < 2
    warning( ...
        'Expected at least two passive-movie repeats in %s; found %d.', ...
        videoTimesPath, size(videoTimes, 2));
    out = struct([]);
    return
end


%% Map all movie repeats onto neural frames

repeatIdx = getAllRepeatFrameIndices(frameTimes, videoTimes);

nRepeats = numel(repeatIdx.start);

if nRepeats < 2
    warning('At least two passive-movie repeats are required.');
    out = struct([]);
    return
end

repeatLengths = repeatIdx.stop - repeatIdx.start + 1;
nFrames = min(repeatLengths);

if max(repeatLengths) - min(repeatLengths) > opts.DurationToleranceFrames
    warning(['Passive-movie repeat durations differ: %s frames. ' ...
        'Truncating all repeats to %d frames.'], ...
        mat2str(repeatLengths), nFrames);
end

repeatFrameIdx = zeros(nFrames, nRepeats);

for iRepeat = 1:nRepeats
    repeatFrameIdx(:,iRepeat) = ...
        repeatIdx.start(iRepeat) : repeatIdx.start(iRepeat) + nFrames - 1;
end

time = frameTimes(repeatFrameIdx(:,1));
time = time - time(1);


%% Process each FOV

out = repmat(struct( ...
    'fov', [], ...
    'fovFolder', '', ...
    'tracePath', '', ...
    'outputPath', '', ...
    'nROIs', [], ...
    'nFrames', nFrames, ...
    'r', [], ...
    'saved', false), 1, numel(fovs));

for iFOV = 1:numel(fovs)

    fovNum = fovs(iFOV);
    fovFolder = fullfile(alfPath, sprintf('FOV_%02d', fovNum));

    tracePath = fullfile(fovFolder, opts.TraceFile);
    outputPath = fullfile(fovFolder, opts.OutputFile);

    out(iFOV).fov = fovNum;
    out(iFOV).fovFolder = fovFolder;
    out(iFOV).tracePath = tracePath;
    out(iFOV).outputPath = outputPath;

    if ~isfile(tracePath)
        warning('Trace file not found for FOV %d: %s', fovNum, tracePath);
        continue
    end

    if isfile(outputPath) && ~opts.overwriteExisting
        fprintf( ...
            'FOV %02d: %s already exists; skipping.\n', ...
            fovNum, opts.OutputFile);

        try
            rExisting = readNPY(outputPath);
            out(iFOV).r = rExisting(:);
            out(iFOV).nROIs = numel(rExisting);
        catch ME
            warning( ...
                'Could not read existing output %s: %s', ...
                outputPath, ME.message);
        end

        continue
    end

    fprintf('FOV %02d: loading traces... ', fovNum);

    traces = double(readNPY(tracePath));
    traces = orientTracesAsTimeByROI( ...
        traces, numel(frameTimes), tracePath);

    maxRequiredFrame = max(repeatFrameIdx(:));

    if size(traces,1) < maxRequiredFrame
        warning( ...
            ['Trace file for FOV %d has only %d frames, but the passive-' ...
            'movie repeat indices require frame %d. Skipping this FOV.'], ...
            fovNum, size(traces,1), maxRequiredFrame);
        continue
    end

    [r, pairwiseR, repeatPairs] = ...
        meanPairwiseRepeatCorrelation(traces, repeatFrameIdx);

    r = r(:);

    out(iFOV).nROIs = size(traces,2);
    out(iFOV).r = r;
    out(iFOV).pairwiseR = pairwiseR;
    out(iFOV).repeatPairs = repeatPairs;
    out(iFOV).nRepeats = nRepeats;

    fprintf('computing %d ROI correlations... ', numel(r));

    writeNPY(r, outputPath);

    out(iFOV).saved = true;

    fprintf('saved %s.\n', outputPath);

    if opts.MakePlots
        plotReliabilityDistribution(r, fovNum, datpath);
    end
end


%% Session-level metadata

for iFOV = 1:numel(out)
    out(iFOV).frameTimesPath = frameTimesPath;
    out(iFOV).videoTimesPath = videoTimesPath;
    out(iFOV).repeatFrameIdx = repeatFrameIdx;
    out(iFOV).repeatBounds = repeatIdx;
    out(iFOV).nRepeats = nRepeats;
    out(iFOV).time = time;
end

end


function frameTimesPath = findFrameTimesFile(alfPath, fovs, fileName)
% Prefer FOV_00, then search the requested FOVs.

frameTimesPath = '';

preferredPath = fullfile(alfPath, 'FOV_00', fileName);

if isfile(preferredPath)
    frameTimesPath = preferredPath;
    return
end

for fovNum = fovs
    candidate = fullfile( ...
        alfPath, sprintf('FOV_%02d', fovNum), fileName);

    if isfile(candidate)
        frameTimesPath = candidate;
        return
    end
end

end


function traces = orientTracesAsTimeByROI(traces, nFrameTimes, tracePath)
% Ensure traces are nFrames x nROIs.

if ndims(traces) ~= 2
    error( ...
        'Expected a 2-D trace matrix in %s; got size %s.', ...
        tracePath, mat2str(size(traces)));
end

if size(traces, 1) == nFrameTimes
    return
end

if size(traces, 2) == nFrameTimes
    warning( ...
        'Transposing traces in %s to obtain nFrames x nROIs.', ...
        tracePath);

    traces = traces';
    return
end

error( ...
    ['Neither dimension of %s matches the number of frame times. ' ...
     'Trace size = %s; number of frame times = %d.'], ...
    tracePath, mat2str(size(traces)), nFrameTimes);

end


function filePath = findFileRecursive(rootDir, fileName)

matches = dir(fullfile(rootDir, '**', fileName));

if isempty(matches)
    filePath = '';
    return
end

if numel(matches) > 1
    matchPaths = fullfile({matches.folder}, {matches.name});

    warning( ...
        ['Found multiple copies of %s. Using the first:\n%s\n' ...
         'All matches:\n%s'], ...
        fileName, matchPaths{1}, strjoin(matchPaths, '\n'));
end

filePath = fullfile(matches(1).folder, matches(1).name);

end


function repeatIdx = getAllRepeatFrameIndices(frameTimes, videoTimes)
%GETALLREPEATFRAMEINDICES Find neural-frame bounds for every movie repeat.

frameTimes = frameTimes(:);
nRepeats = size(videoTimes,2);

repeatIdx.start = nan(1,nRepeats);
repeatIdx.stop  = nan(1,nRepeats);

for iRepeat = 1:nRepeats

    thisVideoTimes = videoTimes(:,iRepeat);
    thisVideoTimes = thisVideoTimes(isfinite(thisVideoTimes));

    if isempty(thisVideoTimes)
        continue
    end

    [~,repeatIdx.start(iRepeat)] = ...
        min(abs(frameTimes - thisVideoTimes(1)));

    [~,repeatIdx.stop(iRepeat)] = ...
        min(abs(frameTimes - thisVideoTimes(end)));
end

valid = isfinite(repeatIdx.start) & isfinite(repeatIdx.stop);
repeatIdx.start = repeatIdx.start(valid);
repeatIdx.stop  = repeatIdx.stop(valid);

if any(repeatIdx.stop < repeatIdx.start)
    error('At least one movie-repeat end precedes its start.');
end

end


function r = columnwisePearson(X, Y)
%COLUMNWISEPEARSON Pearson correlation between corresponding columns.

if ~isequal(size(X), size(Y))
    error( ...
        'Repeat matrices must have identical sizes; got %s and %s.', ...
        mat2str(size(X)), mat2str(size(Y)));
end

% This matches the behavior of ordinary Pearson correlation when there
% are no NaNs. ROIs containing NaNs are returned as NaN.
hasNaN = any(~isfinite(X), 1) | any(~isfinite(Y), 1);

X = X - mean(X, 1);
Y = Y - mean(Y, 1);

numerator = sum(X .* Y, 1);
denominator = sqrt(sum(X.^2, 1) .* sum(Y.^2, 1));

r = numerator ./ denominator;

r(denominator == 0) = NaN;
r(hasNaN) = NaN;

end


function plotReliabilityDistribution(r, fovNum, datpath)

validR = r(isfinite(r));

figure;
histogram( ...
    validR, ...
    -1:0.025:1, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

xlabel('Passive-movie cross-repeat Pearson r');
ylabel('ROI count');
title( ...
    sprintf('%s | FOV %02d', datpath, fovNum), ...
    'Interpreter', 'none');

box off;

end


function [meanR, pairwiseR, repeatPairs] = ...
    meanPairwiseRepeatCorrelation(traces, repeatFrameIdx)
%MEANPAIRWISEREPEATCORRELATION
% Compute every pairwise repeat correlation for every ROI.
%
% Inputs
% ------
% traces
%   nTotalFrames x nROIs
%
% repeatFrameIdx
%   nFramesPerRepeat x nRepeats
%
% Outputs
% -------
% meanR
%   1 x nROIs mean correlation across all repeat pairs
%
% pairwiseR
%   nPairs x nROIs correlation for each repeat pair
%
% repeatPairs
%   nPairs x 2 repeat indices corresponding to rows of pairwiseR

nRepeats = size(repeatFrameIdx,2);
nROIs = size(traces,2);

repeatPairs = nchoosek(1:nRepeats,2);
nPairs = size(repeatPairs,1);

pairwiseR = nan(nPairs,nROIs);

for iPair = 1:nPairs

    iRepeat1 = repeatPairs(iPair,1);
    iRepeat2 = repeatPairs(iPair,2);

    X = traces(repeatFrameIdx(:,iRepeat1),:);
    Y = traces(repeatFrameIdx(:,iRepeat2),:);

    pairwiseR(iPair,:) = columnwisePearson(X,Y);
end

meanR = mean(pairwiseR,1,'omitnan');

% Leave ROIs with no valid pairwise correlations as NaN
meanR(all(isnan(pairwiseR),1)) = NaN;

end