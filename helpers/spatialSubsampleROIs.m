function [Tsub, keepIdx] = spatialSubsampleROIs(T, varargin)
% spatialSubsampleROIs  Pseudo-randomly subsample ROIs for flatter spatial density.
%
% Bins ROIs by ML/AP and keeps up to maxPerBin ROIs per spatial bin.
% Sparse bins keep all ROIs. Dense bins are randomly downsampled.
%
% Usage:
%   Tsub = spatialSubsampleROIs(T);
%   Tsub = spatialSubsampleROIs(T, 'binSize', 100, 'maxPerBin', 50, 'seed', 1);

p = inputParser;
addParameter(p, 'binSize', 100, @(x)isnumeric(x) && isscalar(x));      % microns
addParameter(p, 'maxPerBin', 100, @(x)isnumeric(x) && isscalar(x));     % ROIs per bin
addParameter(p, 'seed', 1, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'xLimits', [], @(x)isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'yLimits', [], @(x)isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'verbose', false, @islogical);
parse(p, varargin{:});
opts = p.Results;

ML = T.ML;
AP = T.AP;

valid = isfinite(ML) & isfinite(AP);

if isempty(opts.xLimits)
    xLimits = [floor(min(ML(valid))/opts.binSize)*opts.binSize, ...
               ceil(max(ML(valid))/opts.binSize)*opts.binSize];
else
    xLimits = opts.xLimits;
end

if isempty(opts.yLimits)
    yLimits = [floor(min(AP(valid))/opts.binSize)*opts.binSize, ...
               ceil(max(AP(valid))/opts.binSize)*opts.binSize];
else
    yLimits = opts.yLimits;
end

xEdges = xLimits(1):opts.binSize:xLimits(2);
yEdges = yLimits(1):opts.binSize:yLimits(2);

if xEdges(end) < xLimits(2)
    xEdges(end+1) = xLimits(2);
end
if yEdges(end) < yLimits(2)
    yEdges(end+1) = yLimits(2);
end

% Bin each ROI
[~,~,xBin] = histcounts(ML, xEdges);
[~,~,yBin] = histcounts(AP, yEdges);

inRange = valid & xBin > 0 & yBin > 0;

keepIdx = false(height(T),1);

rng(opts.seed);

% Convert 2D bin IDs into one group ID
nY = numel(yEdges) - 1;
binID = nan(height(T),1);
binID(inRange) = sub2ind([numel(xEdges)-1, nY], xBin(inRange), yBin(inRange));

%collect some stats
nBinsSubsampled = 0;
nRemoved = 0;
maxOccupancy = 0;

uBins = unique(binID(inRange));

for i = 1:numel(uBins)

    idx = find(binID == uBins(i));
    n = numel(idx);

    maxOccupancy = max(maxOccupancy, n);

    if n <= opts.maxPerBin

        keepIdx(idx) = true;

    else

        nBinsSubsampled = nBinsSubsampled + 1;
        nRemoved = nRemoved + (n - opts.maxPerBin);

        keepIdx(idx(randperm(n, opts.maxPerBin))) = true;

    end

end

Tsub = T(keepIdx,:);

if opts.verbose

    nIn = height(T);
    nOut = height(Tsub);

    if nRemoved == 0

        fprintf(['spatialSubsampleROIs: no subsampling required ' ...
            '(%d ROIs, %d occupied bins, max occupancy %d).\n'], ...
            nIn, numel(uBins), maxOccupancy);

    else

        fprintf(['spatialSubsampleROIs: %d -> %d ROIs ' ...
            '(removed %d, %.1f%%), %d/%d bins subsampled, ' ...
            'max occupancy %d, binSize=%g µm, maxPerBin=%d.\n'], ...
            nIn, ...
            nOut, ...
            nRemoved, ...
            100*nRemoved/nIn, ...
            nBinsSubsampled, ...
            numel(uBins), ...
            maxOccupancy, ...
            opts.binSize, ...
            opts.maxPerBin);

    end

end

end