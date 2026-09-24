function thresholds = get_neuralQMThresholds()
% Default thresholds for neural quality metrics.
%
% Each metric has:
%   value  : threshold value
%   passIf : comparison determining whether the metric passes
%            ('<', '<=', '>', or '>=')

thresholds.noiseLevel            = struct('value',30,   'passIf','<');
thresholds.mean                  = struct ('value',0,   'passIf','>');
thresholds.std                   = struct('value',1,    'passIf','>');
thresholds.skew                  = struct('value',0.1,  'passIf','>');
thresholds.var                   = struct('value',0,    'passIf','>');
thresholds.signalFraction        = struct('value',0.1,  'passIf','>'); % Fraction of detrended neuropil-corrected variance exceeding estimated frame-level noise variance
thresholds.snrTransient          = struct('value',5,    'passIf','>'); % Largest transient amplitude relative to ROI's own frame-level noise SD
thresholds.resNeuropilR2         = struct('value',0.9,  'passIf','<'); % remove ROIs that almost perfectly correlate with neuropil
thresholds.saturationRatio       = struct('value',0.001,'passIf','<'); % remove ROIs whose F is 'clipped'
thresholds.upperRebound          = struct('value',0.2,  'passIf','<'); % remove ROIs that look saturated
thresholds.s2p_footprint         = struct('value',Inf,  'passIf','<');
thresholds.s2p_mrs               = struct('value',Inf,  'passIf','<');
thresholds.s2p_mrs0              = struct('value',1.5,  'passIf','>'); %remove very small ROIs
thresholds.s2p_compact           = struct('value',1,    'passIf','>=');
thresholds.s2p_solidity          = struct('value',0,    'passIf','>');
thresholds.s2p_npix              = struct('value',Inf,  'passIf','<');
thresholds.s2p_npix_soma         = struct('value',10,   'passIf','>'); 
thresholds.s2p_radius            = struct('value',Inf,  'passIf','<');
thresholds.s2p_aspect_ratio      = struct('value',Inf,  'passIf','<');
thresholds.s2p_npix_norm_no_crop = struct('value',Inf,  'passIf','<');
thresholds.s2p_npix_norm         = struct('value',Inf,  'passIf','<');
thresholds.s2p_skew              = struct('value',0,    'passIf','>'); %remove ROIs with negative skew
thresholds.s2p_std               = struct('value',-Inf, 'passIf','>');


end