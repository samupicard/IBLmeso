function bias = computeExponentialBias(stims, tau)
% Compute exponentially weighted average of recent stimuli, 
% with a time constant tau.
%
% Inputs:
%   stims - Vector (floats between 0 and 1, or nan) of size [1, T] or [T, 1].
%   tau   - Time constant for exponential weighting.
%
% Output:
%   bias  - Vector of the same size as 'stims', representing the
%           exponentially weighted average.
%
% Samuel Picard (Dec 2024)

% Ensure sides is a column vector for consistency
stims = stims(:);

% Initialize parameters
T = length(stims);
alpha = 1 - exp(-1 / tau); % Compute alpha
bias = zeros(T, 1); % Preallocate the bias vector
bias(1) = 0.5; % Initial value

% Compute the exponentially weighted average
for t = 2:T
    if isnan(stims(t)) %if stim is 'nan' (i.e. 0-contrast), keep bias the same
        bias(t,:) = bias(t-1,:);
    else
        bias(t) = (1 - alpha) * bias(t-1) + alpha * stims(t);
    end
end

end