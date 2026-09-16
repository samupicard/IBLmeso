function roiIdx = plotQMquantileTraces(traces, metricVals, varargin)
% plotQMquantileTraces
%
% Plot example ROI traces spanning the distribution of a quality metric.
% One ROI is selected nearest each of 10 evenly spaced quantiles.
%
% Inputs
%   traces      : T x N matrix of ROI traces
%   metricVals  : N x 1 vector of quality metric values
%
% Optional name-value pairs
%   'times'     : T x 1 time vector (default = frame number)
%   'mask'      : N x 1 logical mask selecting ROIs (default = all)
%   'nQuantiles': number of example traces (default = 10)
%   'metricName': name used in figure labels (default = 'metric')
%   'traceName' : y-axis label (default = 'trace')
%
% Output
%   roiIdx      : indices of selected ROIs

p = inputParser;
addParameter(p,'times',[],@isnumeric);
addParameter(p,'mask',true(size(traces,2),1),@islogical);
addParameter(p,'nQuantiles',5,@(x)isnumeric(x) && isscalar(x) && x>0);
addParameter(p,'metricName','metric',@(x)ischar(x) || isstring(x));
addParameter(p,'traceName','trace',@(x)ischar(x) || isstring(x));
p.parse(varargin{:});

times = p.Results.times;
mask = p.Results.mask(:);
nQuantiles = p.Results.nQuantiles;
metricName = p.Results.metricName;
traceName = p.Results.traceName;

metricVals = metricVals(:);
nROIs = numel(metricVals);

% Make traces T x N
if size(traces,2) == nROIs
    % already T x N
elseif size(traces,1) == nROIs
    traces = traces';
else
    error(['Neither dimension of traces matches the number of metric values.\n' ...
        'size(traces) = [%d %d], numel(metricVals) = %d'], ...
        size(traces,1),size(traces,2),nROIs);
end

% Time axis
if isempty(times)
    times = (1:size(traces,1))';
    xLabel = 'frame';
else
    times = times(:);

    if numel(times) ~= size(traces,1)
        error(['Time vector length does not match number of frames.\n' ...
            'numel(times) = %d, number of frames = %d'], ...
            numel(times),size(traces,1));
    end

    xLabel = 'time (s)';
end

% Valid candidate ROIs
valid = mask & isfinite(metricVals);

candidateIdx = find(valid);
candidateVals = metricVals(valid);

if numel(candidateIdx) < nQuantiles
    error('Only %d valid ROIs available for %d requested quantiles.', ...
        numel(candidateIdx),nQuantiles);
end

%% choose ROIs nearest evenly spaced quantiles

% Use centres of 10 equal-probability bins:
% 5th, 15th, ..., 95th percentile for nQuantiles = 10
q = ((1:nQuantiles)-0.5) / nQuantiles * 100;
targetVals = prctile(candidateVals,q);

roiIdx = nan(1,nQuantiles);

for i = 1:nQuantiles

    [~,j] = min(abs(candidateVals-targetVals(i)));
    roiIdx(i) = candidateIdx(j);

end

%% plot

figure;
t = tiledlayout(nQuantiles,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for i = 1:nQuantiles

    ax = nexttile;

    plot(times,traces(:,roiIdx(i)));

    ylabel(sprintf('Q%d',i));

    xlim([0,5000]);

    title(sprintf('%s = %.3g   |   ROI %d   |   %.0fth percentile', ...
        metricName,metricVals(roiIdx(i)),roiIdx(i),q(i)), ...
        'Interpreter','none', ...
        'FontWeight','normal');

    box off;

    if i < nQuantiles
        set(ax,'XTickLabel',[]);
    else
        xlabel(xLabel);
    end

end

ylabel(t,traceName,'Interpreter','none');

title(t,sprintf('%s quantile example traces',metricName), ...
    'Interpreter','none');

end