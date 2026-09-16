function thresholds = get_neuralQMThresholds()
% Default thresholds for neural quality metrics

thresholds.noiseLevel      = 100;
thresholds.mean            = 0;
thresholds.std             = 1;
thresholds.skew            = 0;
thresholds.var             = 0;
thresholds.snrVar          = 1;
thresholds.snrTransient    = 5;
thresholds.residualNeuropilR2 = 1;
thresholds.saturationRatio = 0.001;
thresholds.upperRebound    = 0.1;

end