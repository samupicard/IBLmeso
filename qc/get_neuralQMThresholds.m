function thresholds = get_neuralQMThresholds()
% Default thresholds for neural quality metrics.
%
% Each metric has:
%   value  : threshold value
%   passIf : comparison determining whether the metric passes
%            ('<', '<=', '>', or '>=')

thresholds.noiseLevel          = struct('value',30,    'passIf','<');
thresholds.mean                = struct('value',50,    'passIf','>');
thresholds.std                 = struct('value',1,     'passIf','>');
thresholds.skew                = struct('value',0.1,   'passIf','>');
thresholds.var                 = struct('value',0,     'passIf','>');
thresholds.signalFraction      = struct('value',0.2,   'passIf','>'); % Fraction of detrended neuropil-corrected variance exceeding estimated frame-level noise variance
thresholds.snrTransient        = struct('value',3,     'passIf','>'); % Large transient amplitude relative to ROI's own frame-level noise SD
thresholds.resNeuropilR2       = struct('value',0.9,   'passIf','<');
thresholds.saturationRatio     = struct('value',0.001, 'passIf','<');
thresholds.upperRebound        = struct('value',0.2,   'passIf','<');

end