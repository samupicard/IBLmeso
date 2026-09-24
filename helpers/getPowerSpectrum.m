function [powerSpectrum,freq] = getPowerSpectrum(timeseries,Fs,varargin)
%GETPOWERSPECTRUM Multitaper power spectrum of a timecourse.
%
% [powerSpectrum,freq] = getPowerSpectrum(timeseries,Fs)
%
% INPUTS
%   timeseries   : timecourse, vector
%   Fs           : sampling rate of timeseries, in Hz
%
% OPTIONAL NAME-VALUE INPUTS
%   'timeBandwidth' : time-bandwidth product used by PMTM
%                     (default 4)
%
% OUTPUTS
%   powerSpectrum : one-sided power spectral density estimate
%                   (power / Hz)
%   freq          : frequency corresponding to each spectral bin, in Hz
%
% The mean is removed before estimating the spectrum to prevent the DC
% component from dominating the estimate.

p = inputParser;

addParameter(p,'timeBandwidth',6,...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p,varargin{:});
opt = p.Results;


%% Prepare signal

x = double(timeseries(:));

valid = isfinite(x);

if ~all(valid)
    x = x(valid);
end

assert(numel(x) >= 2,...
    'timeseries must contain at least two finite samples.');

assert(isnumeric(Fs) && isscalar(Fs) && Fs > 0,...
    'Fs must be a positive scalar.');


%% Remove DC component

x = x - mean(x);


%% Multitaper PSD

[powerSpectrum,freq] = pmtm( ...
    x, ...
    opt.timeBandwidth, ...
    [], ...
    Fs);

powerSpectrum = powerSpectrum(:);
freq = freq(:);

end