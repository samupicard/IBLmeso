function H = computeStimulusHistory(cDiff,tau,sessionID,mode)
%COMPUTESTIMULUSHISTORY Exponentially weighted stimulus-side history.
%
% H = computeStimulusHistory(cDiff,tau,sessionID,mode)
%
% Inputs:
%   cDiff     : signed contrast difference (right - left)
%   tau       : exponential decay constant, in nonzero-contrast trials
%   sessionID : session identifier for each trial (optional)
%   mode      : 'side' or 'contrast' (default: 'side')
%
% Outputs:
%   H         : Nx1 history regressor, between -1 and 1.
%
% Properties:
% - Only preceding stimuli contribute to each trial.
% - Zero-contrast trials do not update or decay the history.
% - History resets at session boundaries.
% - Before the first informative stimulus, H = 0.

if nargin<3 || isempty(sessionID)
    sessionID = ones(numel(cDiff),1);
end
if nargin<4 || isempty(mode)
    mode = 'side';
end

cDiff = double(cDiff(:));
sessionID = sessionID(:);

assert(tau>0 && isfinite(tau),'tau must be positive.');
assert(numel(sessionID)==numel(cDiff), ...
    'sessionID and cDiff must have the same length.');

alpha = exp(-1/tau);

H = zeros(size(cDiff));

numerator = 0;
denominator = 0;

for t = 1:numel(cDiff)

    % Reset at session boundaries
    if t==1 || ~isequaln(sessionID(t),sessionID(t-1))
        numerator = 0;
        denominator = 0;
    end

    % History preceding the current stimulus
    if denominator>0
        H(t) = numerator/denominator;
    end

    % Update after computing H(t)
    x = cDiff(t);

    if isfinite(x) && x~=0

        switch lower(mode)
            case 'side'
                value = sign(x);
                weight = 1;

            case 'contrast'
                value = sign(x);
                weight = abs(x);

            otherwise
                error('Unknown mode: %s',mode);
        end

        numerator = alpha*numerator + weight*value;
        denominator = alpha*denominator + weight;
    end
end

end
