function thresholds = get_neuralQMThresholds()
% Default thresholds for neural quality metrics.
%
% Each metric has:
%   value  : threshold value
%   passIf : comparison determining whether the metric passes
%            ('lt', 'le', 'gt', or 'ge')

thresholds.noiseLevel          = struct('value',100,   'passIf','<');
thresholds.mean                = struct('value',0,     'passIf','>');
thresholds.std                 = struct('value',1,     'passIf','>=');
thresholds.skew                = struct('value',0,     'passIf','>');
thresholds.var                 = struct('value',0,     'passIf','>');
thresholds.snrVar              = struct('value',2,     'passIf','>');
thresholds.snrTransient        = struct('value',5,     'passIf','>');
thresholds.residualNeuropilR2  = struct('value',1,     'passIf','<=');
thresholds.saturationRatio     = struct('value',0.001, 'passIf','<');
thresholds.upperRebound        = struct('value',0.1,   'passIf','<');
end