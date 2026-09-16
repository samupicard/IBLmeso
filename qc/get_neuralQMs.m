function [neuralQMs, fovQMs] = get_neuralQMs(F,Fneu,varargin)
%GET_NEURALQMS Compute neural quality metrics from raw fluorescence traces.
%
% F and Fneu must be nTimepoints x nROIs, with nROIs < nTimepoints.
%
% INPUTS
% F                 : raw ROI fluorescence traces (nTimepoints x nROIs),
%                     e.g. mpci.ROIActivityF
% Fneu              : raw neuropil fluorescence traces (nTimepoints x nROIs),
%                     e.g. mpci.ROINeuropilActivityF
%
% Optional positional input:
% times             : frame times in seconds (nTimepoints x 1). If supplied,
%                     frame rate Fr is estimated from these times.
%
% Optional name-value parameters:
% 'badframes'       : logical vector marking frames to exclude
%                     (default = all false)
% 'iscell'          : logical vector marking ROIs classified as cells
%                     (default = all true)
% 'Fr'              : frame rate in Hz if 'times' is not supplied
%                     (default = 7)
% 'F0prctile'       : percentile used to compute baseline fluorescence F0
%                     for dF/F; [] disables normalization (default = 20)
% 'transientPrctile': percentile used to estimate transient amplitude
%                     (default = 99.9)
% 'neuropilFactor'  : neuropil subtraction factor (default = 0.7)
% 'saturationTol'   : fraction of each ROI's raw fluorescence range
%                     considered "near maximum" (default = 0.01)
% 'detrendWindow'   : moving-median detrending window in seconds
%                     (default = 60)
%
% OUTPUTS
% neuralQMs         : 1 x nROIs struct with fields:
%
%   noiseLevel          standardized shot-noise level computed from dF/F
%   mean                mean neuropil-corrected fluorescence
%   std                 SD of neuropil-corrected fluorescence
%   skew                skewness of neuropil-corrected fluorescence
%   var                 variance of detrended raw ROI fluorescence
%   snrVar              ratio of detrended ROI fluorescence variance to
%                       detrended neuropil variance
%   snrTransient        transient amplitude divided by robust neuropil noise
%   residualNeuropilR2  squared correlation between detrended,
%                       neuropil-corrected ROI fluorescence and neuropil
%   saturationRatio     fraction of retained frames at or near that ROI's
%                       maximum raw fluorescence
%   upperRebound        largest rise in the coarse raw-fluorescence
%                       histogram after its main peak, normalized by the
%                       main peak height. Values near 0 indicate a roughly
%                       monotonically decreasing upper tail; larger values
%                       indicate a pronounced secondary mode.
%
% fovQMs            : scalar struct containing the mean of each neural QM
%                     across ROIs for which iscell == true.
%
% Samuel Picard (Dec 2024; updated Sep 2026)

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
defaultIscell = true(1,nROIs);
defaultBadframes = false(nTimepoints,1);
defaultSaturationTol = 0.01;
defaultDetrendWindow = 60; % seconds

p = inputParser;

p.addOptional('times', [], ...
    @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x)==nTimepoints))

p.addParameter('badframes', defaultBadframes, ...
    @(x) islogical(x) && isvector(x) && numel(x)==nTimepoints)

p.addParameter('iscell', defaultIscell, ...
    @(x) islogical(x) && isvector(x) && numel(x)==nROIs)

p.addParameter('Fr', defaultFr, ...
    @(x) isnumeric(x) && isscalar(x) && x>0)

p.addParameter('F0prctile', defaultF0prctile, ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=0 && x<=100))

p.addParameter('transientPrctile', defaultTransientPrctile, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=100)

p.addParameter('neuropilFactor', defaultNeuropilFactor, ...
    @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1)

p.addParameter('saturationTol', defaultSaturationTol, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1)

p.addParameter('detrendWindow', defaultDetrendWindow, ...
    @(x) isnumeric(x) && isscalar(x) && x>0)

p.parse(varargin{:})

times = p.Results.times;
badframes = p.Results.badframes(:);
iscell = p.Results.iscell(:)';
Fr = p.Results.Fr;
F0prctile = p.Results.F0prctile;
transientPrctile = p.Results.transientPrctile;
neuropilFactor = p.Results.neuropilFactor;
saturationTol = p.Results.saturationTol;
detrendWindow = p.Results.detrendWindow;

%% Re-work some arguments

if isempty(times)

    if Fr == defaultFr
        warning('Neither frame rate nor frame times provided, assuming %.1f Hz.',Fr)
    end

else

    times = times(:);
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
detrendWindowFrames = max(3,round(detrendWindow*Fr));

F_baseline = movmedian(F,detrendWindowFrames,1, ...
    'omitmissing','Endpoints','shrink');

Fneu_baseline = movmedian(Fneu,detrendWindowFrames,1, ...
    'omitmissing','Endpoints','shrink');

F_detrended = F - F_baseline;
Fneu_detrended = Fneu - Fneu_baseline;

%% Saturation QC

maxF = max(F,[],1);
minF = min(F,[],1);
F_range = maxF - minF;

% Per-ROI threshold defining "at or near maximum"
saturationThreshold = maxF - saturationTol.*F_range;

% Fraction of retained frames at/above that threshold
saturationRatios = mean(F >= saturationThreshold,1);

% Undefined for ROIs with no fluorescence range
badRange = F_range <= 0;
saturationRatios(badRange) = NaN;

%% Compute other QMs

% Standardized shot noise level
noiseLevels = median(abs(diff(dFF,1,1)),1,'omitnan') / sqrt(Fr);

% Mean neuropil-corrected fluorescence
means = mean(F_npc,1);

% Standard deviation of neuropil-corrected fluorescence
stds = std(F_npc,0,1);

% Skewness of neuropil-corrected fluorescence
skews = skewness(F_npc,0,1);

% Variance of detrended ROI fluorescence
vars = var(F_detrended,1,1);

% Ratio of detrended ROI fluorescence variance to neuropil variance
snrs = vars ./ var(Fneu_detrended,1,1);

% Transient amplitude and transient SNR relative to neuropil fluctuations
F_npc_detrended = F_detrended - neuropilFactor*Fneu_detrended;

transientAmp = ...
    prctile(F_npc_detrended,transientPrctile,1) - ...
    median(F_npc_detrended,1);

neuropilNoise = 1.4826 * mad(Fneu_detrended,1,1);
% 1.4826 scales MAD to approximately SD for a Gaussian distribution

transientSNRs = transientAmp ./ neuropilNoise;

% Squared correlation between detrended neuropil-corrected F and neuropil
F_centered = F_npc_detrended - mean(F_npc_detrended,1);
Fneu_centered = Fneu_detrended - mean(Fneu_detrended,1);

r = sum(F_centered.*Fneu_centered,1) ./ ...
    sqrt(sum(F_centered.^2,1).*sum(Fneu_centered.^2,1));

residualNeuropilR2 = r.^2;

%% Compute upper rebound metric (measure of possible saturation)

nBins = 50;
lowPrctile = 1;
highPrctile = 99.9;

% Robustly normalize raw ROI fluorescence
lo = prctile(F,lowPrctile,1);
hi = prctile(F,highPrctile,1);

F_norm = (F-lo) ./ (hi-lo);

% Histogram all ROIs using common normalized bins
edges = linspace(0,1,nBins+1);
roiIdx = repmat(1:nROIs,size(F,1),1);

counts = histcounts2( ...
    F_norm(:),roiIdx(:), ...
    edges,0.5:1:(nROIs+0.5));

% Normalize histogram counts within each ROI
counts = counts ./ sum(counts,1);

% Find main histogram peak
[peakHeight,peakIdx] = max(counts,[],1);

% Ignore bins before the main peak
binIdx = (1:nBins)';
tailMask = binIdx >= peakIdx;

tailCounts = counts;
tailCounts(~tailMask) = NaN;

% Largest subsequent rise above the running minimum
runningMin = cummin(tailCounts,1,'omitnan');
rebound = tailCounts - runningMin;

upperRebound = ...
    max(rebound,[],1,'omitnan') ./ (peakHeight + eps);

% Undefined for ROIs with negligible robust fluorescence range
badRange = (hi-lo) <= 0;
upperRebound(badRange) = NaN;

%% Return neural quality metrics

neuralQMs = struct([]);

for iROI = 1:nROIs

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

fovQMs = struct( ...
    'noiseLevel',         mean(noiseLevels(iscell),'omitnan'), ...
    'mean',               mean(means(iscell),'omitnan'), ...
    'std',                mean(stds(iscell),'omitnan'), ...
    'skew',               mean(skews(iscell),'omitnan'), ...
    'var',                mean(vars(iscell),'omitnan'), ...
    'snrVar',             mean(snrs(iscell),'omitnan'), ...
    'snrTransient',       mean(transientSNRs(iscell),'omitnan'), ...
    'residualNeuropilR2', mean(residualNeuropilR2(iscell),'omitnan'), ...
    'saturationRatio',    mean(saturationRatios(iscell),'omitnan'), ...
    'upperRebound',       mean(upperRebound(iscell),'omitnan'));

end