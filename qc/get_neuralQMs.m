function [neuralQMs, fovQMs] = get_neuralQMs(F,Fneu,varargin)
%compute neural quality metrics based on raw fluorescence traces
%
% F is raw fluorescence trace (nTimepoints by nROIs), e.g. mpci.ROINeuropilActivityF
% Fneu is raw neuropil trace (nTimepoints by nROIs), e.g. mpci.ROIActivityF
% times (optional) is frame times in seconds (nTimepoints by 1), e.g. mpci.times
%
% optional name-value pair parameters:
% 'badframes'       : indices of frames that should be excluded (default = all false)
% 'iscell'          : boolean array with true for cells, and false for not cells (default = all true)
% 'Fr'              : frame rate in Hz, in case it cannot be computed from 'times' (default = 7)
% 'F0prctile'       : percentile to be used for computing baseline fluorescence F0 (default = 20)
% 'TransientPrctile': percentile to be used for computing transient amplitude (default = 99.9)
% 'neuropilFactor'  : factor to multiply neuropil with to get neuropil-corrected trace (default = 0.7)
% 'saturationTol'   : fraction of each ROI's raw fluorescence range considered
%                     "near maximum" (default = 0.01, i.e. top 1%)
% 'detrendWindow'   : window in seconds for moving median detrending (default = 60)
%
% returns
% neuralQMs (1xnROIs struct) with the following fields:
%   noiseLevel      : standardized shot noise level
%   mean            : time-averaged raw fluorescence
%   std             : standard deviation of neuropil-corrected activity
%   skew            : skewness of neuropil-corrected activity
%   saturationRatio : fraction of all recording frames for which raw
%                     fluorescence is at or near that ROI's maximum.
%   upperRebound    : how strongly the fluorescence histogram rises again 
%                     after its main peak. Values near 0 indicate a roughly 
%                     monotonically decreasing upper tail, while larger 
%                     values indicate a pronounced secondary mode, which 
%                     may reflect saturation-like behaviour.
% 
% fovQMs (1x1 struct) with the same fields, averaged across all ROIs
% for which iscell = true
%
% Samuel Picard (Dec 2024)

% Validate F and Fneu dimensions
assert(ismatrix(F) && ismatrix(Fneu), ...
    'F and Fneu must both be 2D matrices.');

assert(isequal(size(F), size(Fneu)), ...
    'F and Fneu must have identical dimensions.');

[nTimepoints, nROIs] = size(F);

assert(nROIs < nTimepoints, ...
    ['F and Fneu must be nTimepoints x nROIs, with nROIs < nTimepoints. ' ...
     'Received size %d x %d.'], nTimepoints, nROIs);

defaultF0prctile = 20;
defaultTransientPrctile = 99.9;
defaultNeuropilFactor = 0.7;
defaultFr = 7;
defaultIscell = true(1,size(F,2));
defaultBadframes = false(1,size(F,1));
defaultSaturationTol = 0.01;
defaultDetrendWindow = 60; % seconds

p = inputParser;

p.addOptional('times', [], @(x) isnumeric(x))
p.addParameter('badframes', defaultBadframes, @islogical)
p.addParameter('iscell', defaultIscell, @islogical)
p.addParameter('Fr', defaultFr, @(x) isnumeric(x) && (x>0))
p.addParameter('F0prctile', defaultF0prctile, @(x) isnumeric(x))
p.addParameter('transientPrctile', defaultTransientPrctile, @(x) isnumeric(x))
p.addParameter('neuropilFactor', defaultNeuropilFactor, ...
    @(x) isnumeric(x) && (x>0) && (x<=1))
p.addParameter('saturationTol', defaultSaturationTol, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1)
p.addParameter('detrendWindow', defaultDetrendWindow, ...
    @(x) isnumeric(x) && isscalar(x) && x>0)

p.parse(varargin{:})

times = p.Results.times;
badframes = p.Results.badframes;
iscell = p.Results.iscell;
Fr = p.Results.Fr;
F0prctile = p.Results.F0prctile;
transientPrctile = p.Results.transientPrctile;
neuropilFactor = p.Results.neuropilFactor;
saturationTol = p.Results.saturationTol;
detrendWindow = p.Results.detrendWindow;

%% Re-work some arguments

if isempty(times)

    if Fr==7
        warning('neither frame rate nor array of frame times provided, assuming 7Hz')
    end

else

    if size(times,1) == 1
        times = times';
    end

    Fr = 1/median(diff(times));

end

%% Keep only good frames

F = F(~badframes,:);
Fneu = Fneu(~badframes,:);

%% Get neuropil-corrected trace and delta F / F0

F_npc = F - neuropilFactor*Fneu;

if ~isempty(F0prctile)
    F0 = prctile(F,F0prctile,1);
    dFF = (F-F0) ./ F0 * 100;
else
    dFF = F;
end

%% Detrend raw fluorescence traces

% Remove slow fluorescence drift/bleaching using a moving median.
% The window should be much longer than GCaMP6s calcium transients.
detrendWindowFrames = max(3, round(detrendWindow * Fr));

F_baseline = movmedian(F, detrendWindowFrames, 1, ...
    'omitmissing', 'Endpoints', 'shrink');

Fneu_baseline = movmedian(Fneu, detrendWindowFrames, 1, ...
    'omitmissing', 'Endpoints', 'shrink');

F_detrended = F - F_baseline;
Fneu_detrended = Fneu - Fneu_baseline;

%% Saturation QC

maxF = max(F,[],1);
minF = min(F,[],1);

F_range = maxF - minF;

% Per-ROI threshold defining "at or near maximum"
saturationThreshold = maxF - saturationTol .* F_range;

% Fraction of full recording at/above that threshold
saturationRatios = mean(F >= saturationThreshold, 1);

% take out metric for near 0-range ROIs
badRange = F_range <= 0;
saturationRatios(badRange') = NaN;

%% Compute other QMs

% standardized shot noise level
noiseLevels = nanmedian(abs(diff(dFF,1,1)),1)/sqrt(Fr);

% time-averaged neuropil-corrected activity
means = mean(F_npc,1);

% standard deviation of neuropil-corrected activity
stds = std(F_npc,0,1);

% skewness of neuropil-corrected activity
skews = skewness(F_npc,0,1);

% variance of detrended ROI fluorescence
vars = var(F_detrended,1); 

% snr: variance of detrended ROI fluorescence / variance of ROI neuropil
snrs = vars ./ var(Fneu_detrended,1);

% transientAmp and transientSNR, relative to neuropil fluctuations
F_npc_detrended = F_detrended - neuropilFactor * Fneu_detrended;
transientAmp = prctile(F_npc_detrended, transientPrctile, 1) - median(F_npc_detrended, 1);
neuropilNoise = 1.4826 * mad(Fneu_detrended, 1, 1); %1.4826 scales MAD to be comparable to SD for gaussian distribution
transientSNRs = transientAmp ./ neuropilNoise;

% correlation between detrended centered F_npc and Fneu
F_centered = F_npc_detrended - mean(F_npc_detrended,1);
Fneu_centered = Fneu_detrended - mean(Fneu_detrended,1);
r = sum(F_centered .* Fneu_centered, 1) ./ ...
    sqrt(sum(F_centered.^2,1) .* sum(Fneu_centered.^2,1));
residualNeuropilR2 = r.^2;

%% Compute upper rebound metric (measure of possible saturation)

% Parameters
nBins = 50;
lowPrctile = 1;
highPrctile = 99.9;

% Robustly normalize each ROI to approximately [0,1]
lo = prctile(F_npc, lowPrctile, 1);
hi = prctile(F_npc, highPrctile, 1);

F_norm = (F_npc - lo) ./ (hi - lo);

% Histogram all ROIs using common normalized bins
edges = linspace(0, 1, nBins+1);

nTimepoints = size(F_npc,1);
nROIs = size(F_npc,2);

roiIdx = repmat(1:nROIs, nTimepoints, 1);

counts = histcounts2( ...
    F_norm(:), roiIdx(:), ...
    edges, 0.5:1:(nROIs+0.5));

% counts is nBins x nROIs
counts = counts ./ sum(counts,1);

% Find the main histogram peak for each ROI
[peakHeight, peakIdx] = max(counts, [], 1);

% Ignore bins before the main peak
binIdx = (1:nBins)';

tailMask = binIdx >= peakIdx;

tailCounts = counts;
tailCounts(~tailMask) = NaN;

% Running minimum after the main peak
runningMin = cummin(tailCounts, 1, 'omitnan');

% Amount by which the histogram subsequently rebounds
rebound = tailCounts - runningMin;

% Largest rebound, normalized by the main peak height
upperRebound = max(rebound, [], 1, 'omitnan') ./ (peakHeight + eps);

% Return as nROIs x 1
upperRebound = upperRebound';

% ROIs with essentially no fluorescence range are undefined
badRange = (hi - lo) <= 0;
upperRebound(badRange') = NaN;

%% Return neural quality metrics

neuralQMs = struct([]);

for iROI = 1:size(F,2)

    neuralQMs(iROI).noiseLevel = noiseLevels(iROI);
    neuralQMs(iROI).mean = means(iROI);
    neuralQMs(iROI).std = stds(iROI);
    neuralQMs(iROI).skew = skews(iROI);
    neuralQMs(iROI).var = vars(iROI);
    neuralQMs(iROI).snrVar = snrs(iROI);
    neuralQMs(iROI).snrTransient = transientSNRs(iROI);
    neuralQMs(iROI).residualNeuropilR2 = residualNeuropilR2(iROI);
    neuralQMs(iROI).saturationRatio = saturationRatios(iROI);
    neuralQMs(iROI).upperRebound = upperRebound(iROI);
end

%% Return FOV-wide mean quality metrics

fovQMs = struct(...
    'noiseLevel',        nanmean(noiseLevels(iscell)),...
    'mean',              nanmean(means(iscell)),...
    'std',               nanmean(stds(iscell)),...
    'skew',              nanmean(skews(iscell)),...
    'var',               nanmean(vars(iscell)),...
    'snrVar',            nanmean(snrs(iscell)),...
    'snrTransient',      nanmean(transientSNRs(iscell)),...
    'residualNeuropilR2',nanmean(residualNeuropilR2(iscell)),...
    'saturationRatio',   nanmean(saturationRatios(iscell)),...
    'upperRebound',      nanmean(upperRebound(iscell)));

end