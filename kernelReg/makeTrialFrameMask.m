function inTrialMask = makeTrialFrameMask( ...
    trialsT, frameTimes, preGoCue, postFeedback)
%MAKETRIALFRAMEMASK Identify frames within task-trial evaluation periods.
%
% A frame is included if it lies between:
%
%   goCue_times - preGoCue
%
% and:
%
%   feedback_times + postFeedback
%
% for at least one trial.

arguments
    trialsT table
    frameTimes (:,1) double
    preGoCue (1,1) double {mustBeNonnegative} = 0.6
    postFeedback (1,1) double {mustBeNonnegative} = 1.2
end

frameTimes = frameTimes(:);

trialStart = double(trialsT.goCue_times(:)) - preGoCue;
trialEnd = double(trialsT.feedback_times(:)) + postFeedback;

validTrials = ...
    isfinite(trialStart) & ...
    isfinite(trialEnd) & ...
    trialEnd >= trialStart;

trialStart = trialStart(validTrials);
trialEnd = trialEnd(validTrials);

inTrialMask = false(size(frameTimes));

for trialIndex = 1:numel(trialStart)
    inTrialMask = inTrialMask | ...
        (frameTimes >= trialStart(trialIndex) & ...
         frameTimes <= trialEnd(trialIndex));
end

end