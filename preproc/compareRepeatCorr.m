function results = compareRepeatCorr(sPaths, varargin)
%COMPAREREPEATCORR Cross-repeat neural response reliability for video stimuli.
%
% results = compareRepeatCorr(sPaths, 'MakePlots', true, 'FOVs', [0 1]);
%
% Either suite2p (ALF/FOV_XX) or maskNMF (suite2p/planeX/masknmf_output)
% may be present independently in each FOV. At least 2 video repeats required.
%
% For K repeats, compute the temporal Pearson r for each of the K*(K-1)/2
% pairs. The summary r is tanh(mean(atanh(pairwise r))), omitting NaNs.
% For K=2, this is simply the original cross-repeat Pearson r.
%
% Output for each session:
%   results(i).s2p / .masknmf:
%       .r              1 x nROIs summary reliability
%       .pairR          nPairs x nROIs individual repeat-pair r values
%       .fovs           1 x nROIs originating FOV number
%       .pos            nROIs x 2 position, NaNs if unavailable
%       .traces.repeat1, .repeat2, ... (nFrames x nROIs)
%       .raw.repeat1, .repeat2, ... (neuropil-corrected for suite2p)
%   Missing methods have empty arrays and empty trace/raw structs.
%   .pairIndices is nPairs x 2; .repeatFrameIdx is K x 2 [start end].
%
% All repeat windows are truncated to the shortest repeat, matching the
% original behaviour. This assumes approximately equal frame rates and
% stimulus timing across repeats, without time-warping/interpolation.

p = inputParser;
p.addParameter('TraceFile', 'mpci.ROIActivityDeconvolved.npy', @ischar);
p.addParameter('RawFile', 'mpci.ROIActivityF.npy', @ischar);
p.addParameter('NeuropilFile', 'mpci.ROINeuropilActivityF.npy', @ischar);
p.addParameter('NeuropilWeight', 0.7, @isscalar);
p.addParameter('ClassifierFile', 'mpciROIs.cellClassifier.npy', @ischar);
p.addParameter('PosFile', 'mpciROIs.stackPos.npy', @ischar);
p.addParameter('FrameTimesFile', 'mpci.times.npy', @ischar);
p.addParameter('VideoTimesFile', '_sp_video.times.npy', @ischar);
p.addParameter('DurationToleranceFrames', 1, @isscalar);
p.addParameter('MakePlots', true, @islogical);
p.addParameter('FOVs', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('ClassifierThreshold', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x)));
p.parse(varargin{:});
opts = p.Results;

if ischar(sPaths) || isstring(sPaths)
    sPaths = cellstr(sPaths);
end
results = struct([]);

for iPath = 1:numel(sPaths)
    sPath = sPaths{iPath};
    fprintf('Processing %s\n', sPath);

    [fovS2p, fovMask] = findAvailableFOVs(sPath, opts.TraceFile);
    if ~isempty(opts.FOVs)
        fovS2p = intersect(opts.FOVs, fovS2p, 'stable');
        fovMask = intersect(opts.FOVs, fovMask, 'stable');
        missing = setdiff(opts.FOVs, union(fovS2p, fovMask));
        if ~isempty(missing)
            warning('Requested FOVs missing both methods in %s: %s', ...
                sPath, mat2str(missing));
        end
    end
    fovs = union(fovS2p, fovMask, 'stable');
    if isempty(fovs)
        warning('No usable FOVs in %s. Skipping.', sPath);
        continue
    end

    videoTimesPath = findFileRecursive(fullfile(sPath, 'alf'), opts.VideoTimesFile);
    if isempty(videoTimesPath)
        warning('Could not find %s in %s. Skipping.', opts.VideoTimesFile, sPath);
        continue
    end
    videoTimes = readNPY(videoTimesPath);
    if size(videoTimes, 2) < 2
        warning('Fewer than two video repeats in %s. Skipping.', sPath);
        continue
    end
    nRepeats = size(videoTimes, 2);

    % All FOVs are assumed to have the same frame index/time base.
    frameTimesPath = '';
    for iFOV = 1:length(fovs)
        fov = fovs(iFOV);
        candidate = fullfile(sPath, 'alf', sprintf('FOV_%02d', fov), opts.FrameTimesFile);
        if isfile(candidate)
            frameTimesPath = candidate;
            break
        end
    end
    if isempty(frameTimesPath)
        warning('Frame times missing in %s. Skipping.', sPath);
        continue
    end
    frameTimes = readNPY(frameTimesPath);
    frameTimes = frameTimes(:);
    repIdx = getRepeatFrameIndices(frameTimes, videoTimes);
    lengths = repIdx(:,2) - repIdx(:,1) + 1;
    if any(lengths <= 2)
        warning('Repeat window too short or invalid in %s. Skipping.', sPath);
        continue
    end
    if max(lengths) - min(lengths) > opts.DurationToleranceFrames
        warning('Repeat lengths differ in %s (%s frames); truncating.', ...
            sPath, num2str(lengths(:)'));
    end
    nFrames = min(lengths);
    frameIndices = cell(1, nRepeats);
    for k = 1:nRepeats
        frameIndices{k} = repIdx(k,1):(repIdx(k,1) + nFrames - 1);
    end
    timeAxis = frameTimes(frameIndices{1});
    timeAxis = timeAxis - timeAxis(1);

    % Load methods independently; where both are available for an FOV,
    % retain original count-matching strategy (top suite2p classifier ROIs).
    data = struct('s2p', emptyMethod(), 'masknmf', emptyMethod());
    for iFOV =  1:length(fovs)
        fov = fovs(iFOV);
        s2pPath = fullfile(sPath, 'alf', sprintf('FOV_%02d', fov));
        maskPath = fullfile(sPath, 'suite2p', sprintf('plane%d', fov), 'masknmf_output');
        hasS2p = ismember(fov, fovS2p);
        hasMask = ismember(fov, fovMask);

        if hasMask
            m = loadMethodFOV(maskPath, 'masknmf', opts);
        end
        if hasS2p
            s = loadMethodFOV(s2pPath, 's2p', opts);
            useThreshold = ~isempty(opts.ClassifierThreshold);
            useCountMatch = ~useThreshold && hasMask && size(m.tr,2) < size(s.tr,2);

            if useThreshold || useCountMatch
                classifierPath = fullfile(s2pPath, opts.ClassifierFile);

                if isfile(classifierPath)
                    scores = readNPY(classifierPath);

                    if numel(scores) == size(s.tr,2)
                        if useThreshold
                            keep = find(scores(:) >= opts.ClassifierThreshold);
                        else
                            [~, keep] = maxk(scores(:), size(m.tr,2));
                        end

                        s.tr = s.tr(:,keep);
                        if ~isempty(s.raw), s.raw = s.raw(:,keep); end
                        s.pos = s.pos(keep,:);
                    elseif useThreshold
                        warning('Classifier size mismatch; skipping suite2p FOV.');
                        hasS2p = false;
                    end
                elseif useThreshold
                    warning('Classifier missing; skipping suite2p FOV.');
                    hasS2p = false;
                end
            end
        end
        if hasMask, data.masknmf = appendMethod(data.masknmf, m, fov); end
        if hasS2p, data.s2p = appendMethod(data.s2p, s, fov); end
    end

    pairIndices = nchoosek(1:nRepeats, 2);
    for name = {'s2p', 'masknmf'}
        method = name{1};
        d = data.(method);
        if isempty(d.tr)
            d.r = [];
            d.pairR = [];
            d.traces = struct();
            d.raw = struct();
        else
            if size(d.tr,1) < max(cellfun(@max, frameIndices))
                error('%s traces too short for video repeat frames in %s.', method, sPath);
            end
            d.traces = struct();
            d.raw = struct();
            for k = 1:nRepeats
                field = sprintf('repeat%d', k);
                d.traces.(field) = d.tr(frameIndices{k},:);
                if ~isempty(d.rawSignal)
                    d.raw.(field) = d.rawSignal(frameIndices{k},:);
                end
            end
            d.pairR = nan(size(pairIndices,1), size(d.tr,2));
            for j = 1:size(pairIndices,1)
                a = sprintf('repeat%d', pairIndices(j,1));
                b = sprintf('repeat%d', pairIndices(j,2));
                d.pairR(j,:) = columnwisePearson(d.traces.(a), d.traces.(b));
            end
            d.r = fisherMean(d.pairR);
        end
        d = rmfield(d, {'tr', 'rawSignal'});
        data.(method) = d;
    end

    out.sPath = sPath;
    out.fovs = fovs;
    out.videoTimesPath = videoTimesPath;
    out.nRepeats = nRepeats;
    out.nFrames = nFrames;
    out.repeatFrameIdx = repIdx;
    out.pairIndices = pairIndices;
    out.time = timeAxis;
    out.masknmf = data.masknmf;
    out.s2p = data.s2p;
    results = [results, out]; %#ok<AGROW>

    if opts.MakePlots
        figure; hold on
        legendEntries = {};
        edges = -.1:0.025:1;
        if ~isempty(out.s2p.r)
            histogram(out.s2p.r, edges, 'DisplayStyle', 'stairs', 'LineWidth', 2);
            legendEntries{end+1} = 'suite2p'; %#ok<AGROW>
        end
        if ~isempty(out.masknmf.r)
            histogram(out.masknmf.r, edges, 'DisplayStyle', 'stairs', 'LineWidth', 2);
            legendEntries{end+1} = 'maskNMF'; %#ok<AGROW>
        end
        legend(legendEntries, 'Location', 'best');
        xlabel('Cross-repeat Pearson r (Fisher-z mean across pairs)');
        ylabel('ROI count');
        title(sprintf('%s (%d repeats)', sPath, nRepeats), 'Interpreter', 'none');
    end
end
end

function d = emptyMethod()
d = struct('tr', [], 'rawSignal', [], 'fovs', [], 'pos', [], ...
    'r', [], 'pairR', [], 'traces', struct(), 'raw', struct());
end

function [s2pFOVs, maskFOVs] = findAvailableFOVs(sPath, traceFn)
s2pFOVs = [];
maskFOVs = [];
alfDir = fullfile(sPath, 'alf');
folders = dir(fullfile(alfDir, 'FOV_*'));
for k = 1:numel(folders)
    if ~folders(k).isdir, continue; end
    n = sscanf(folders(k).name, 'FOV_%d');
    if isempty(n), continue; end
    s2pPath = fullfile(alfDir, folders(k).name, traceFn);
    maskPath = fullfile(sPath, 'suite2p', sprintf('plane%d',n), 'masknmf_output', traceFn);
    if isfile(s2pPath), s2pFOVs(end+1) = n; end %#ok<AGROW>
    if isfile(maskPath), maskFOVs(end+1) = n; end %#ok<AGROW>
end
% Mask-only FOVs may not have an ALF/FOV_XX directory.
maskDirs = dir(fullfile(sPath, 'suite2p', 'plane*'));
for k = 1:numel(maskDirs)
    if ~maskDirs(k).isdir, continue; end
    n = sscanf(maskDirs(k).name, 'plane%d');
    if isempty(n), continue; end
    if isfile(fullfile(maskDirs(k).folder, maskDirs(k).name, 'masknmf_output', traceFn))
        maskFOVs(end+1) = n; %#ok<AGROW>
    end
end
s2pFOVs = unique(s2pFOVs, 'stable');
maskFOVs = unique(maskFOVs, 'stable');
end

function d = loadMethodFOV(folder, method, opts)
d.tr = readNPY(fullfile(folder, opts.TraceFile));
d.raw = [];
if isfile(fullfile(folder, opts.RawFile))
    d.raw = readNPY(fullfile(folder, opts.RawFile));
    if strcmp(method, 's2p') && isfile(fullfile(folder, opts.NeuropilFile))
        d.raw = d.raw - opts.NeuropilWeight * readNPY(fullfile(folder, opts.NeuropilFile));
    end
end
nROIs = size(d.tr,2);
d.pos = nan(nROIs, 2);
if isfile(fullfile(folder, opts.PosFile))
    pos = readNPY(fullfile(folder, opts.PosFile));
    if size(pos,1) == nROIs && size(pos,2) >= 2
        d.pos = pos(:,1:2);
    else
        warning('Position metadata size mismatch in %s.', folder);
    end
end
end

function accumulated = appendMethod(accumulated, next, fov)
if ~isempty(accumulated.tr) && size(accumulated.tr,1) ~= size(next.tr,1)
    error('Cannot concatenate FOV %d: trace lengths differ.', fov);
end
oldN = size(accumulated.tr,2);
newN = size(next.tr,2);
accumulated.tr = [accumulated.tr, next.tr];
% If raw is unavailable for any included FOV, keep NaN placeholders.
if ~isempty(accumulated.rawSignal) || ~isempty(next.raw)
    if isempty(accumulated.rawSignal)
        accumulated.rawSignal = nan(size(next.tr,1), oldN);
    end
    if isempty(next.raw)
        next.raw = nan(size(next.tr,1), newN);
    end
    accumulated.rawSignal = [accumulated.rawSignal, next.raw];
end
accumulated.pos = [accumulated.pos; next.pos];
accumulated.fovs = [accumulated.fovs, repmat(fov, 1, newN)];
end

function idx = getRepeatFrameIndices(frameTimes, videoTimes)
% Columns of videoTimes are repeats; first and last rows delimit each.
nRep = size(videoTimes,2);
idx = zeros(nRep,2);
for k = 1:nRep
    [~,idx(k,1)] = min(abs(frameTimes - videoTimes(1,k)));
    [~,idx(k,2)] = min(abs(frameTimes - videoTimes(end,k)));
end
end

function r = columnwisePearson(A, B)
% Columnwise Pearson r using paired finite samples (no toolboxes needed).
valid = isfinite(A) & isfinite(B);
n = sum(valid, 1);
A(~valid) = 0;
B(~valid) = 0;
meanA = sum(A,1) ./ max(n,1);
meanB = sum(B,1) ./ max(n,1);
A = (A - meanA) .* valid;
B = (B - meanB) .* valid;
r = sum(A.*B,1) ./ sqrt(sum(A.^2,1).*sum(B.^2,1));
r(n < 3 | ~isfinite(r)) = NaN;
r = max(-1, min(1, r));
end

function r = fisherMean(pairR)
% Fisher z average, omitting invalid repeat-pair estimates per ROI.
valid = isfinite(pairR);
z = atanh(max(-1+1e-7, min(1-1e-7, pairR)));
z(~valid) = 0;
n = sum(valid,1);
r = tanh(sum(z,1) ./ max(n,1));
r(n == 0) = NaN;
end

function filePath = findFileRecursive(rootDir, fileName)
matches = dir(fullfile(rootDir, '**', fileName));
if isempty(matches)
    filePath = '';
else
    if numel(matches) > 1
        warning('Multiple %s files found; using %s.', fileName, ...
            fullfile(matches(1).folder, matches(1).name));
    end
    filePath = fullfile(matches(1).folder, matches(1).name);
end
end
