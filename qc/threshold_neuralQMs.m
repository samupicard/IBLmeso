function passTable = threshold_neuralQMs(neuralQMs, thresholds)
% Threshold neural quality metrics.
%
% INPUTS
% neuralQMs  : 1 x nROIs struct containing neural quality metrics
% thresholds : optional threshold struct. If omitted, uses defaults from
%              get_neuralQMThresholds
%
% OUTPUT
% passTable  : nROIs x nMetrics logical table; true = metric passed

if nargin < 2 || isempty(thresholds)
    thresholds = get_neuralQMThresholds();
end

qmTable = struct2table(neuralQMs);
metricNames = fieldnames(thresholds);

nROIs = height(qmTable);
nMetrics = numel(metricNames);

pass = false(nROIs,nMetrics);

for iMetric = 1:nMetrics

    metricName = metricNames{iMetric};

    assert(ismember(metricName,qmTable.Properties.VariableNames), ...
        'Metric "%s" is defined in thresholds but missing from neuralQMs.', ...
        metricName);

    vals = qmTable.(metricName);
    thr = thresholds.(metricName).value;

    switch thresholds.(metricName).passIf
        case '<'
            pass(:,iMetric) = vals < thr;
        case '<='
            pass(:,iMetric) = vals <= thr;
        case '>'
            pass(:,iMetric) = vals > thr;
        case '>='
            pass(:,iMetric) = vals >= thr;
        otherwise
            error('Unknown passIf "%s" for metric "%s".', ...
                thresholds.(metricName).passIf,metricName);
    end

end

passTable = array2table(pass,'VariableNames',metricNames);

end