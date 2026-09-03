function [X, predictorInfo, predictors] = buildIBLDesignMatrix( ...
    trialsT, frameTimes, wheelPosition, wheelTimestamps, options)
%BUILDIBLDESIGNMATRIX Construct an IBL task design matrix.
%
% [X, predictorInfo, predictors] = buildIBLDesignMatrix( ...
%     trialsT, frameTimes, wheelPosition, wheelTimestamps)
%
% constructs a design matrix sampled at the two-photon imaging frame times.
%
% DISCRETE EVENT KERNELS
% ----------------------
%
% Each discrete event is expanded over a lag window.
%
%   stimulusLeft:
%       contrastLeft at goCue_times (post)
%
%   stimulusRight:
%       contrastRight at goCue_times (post)
%
%   goCue:
%       unit impulse at goCue_times
%
%   choice:
%       trials.choice at choiceMovement_times (pre)
%
%   action:
%       unit impulse at choiceMovement_times (post)
%
%   feedback:
%       trials.feedbackType at feedback_times
%
% By default:
%
%   goCue          -.6 to .6 s relative to goCue
%   stimulusLeft   0 to 1.2 s relative to goCue
%   stimulusRight  0 to 1.2 s relative to goCue
%   choice         -1 to 0.2 s relative to choiceMovement
%   action       -0.2 to 1.0 s relative to choiceMovement
%   feedback        0 to 1.2 s relative to feedback
%
% EPOCH REGRESSOR
% ---------------
%
%   block, included only if the frame-aligned regressor has three values:
%       probabilityLeft from quiescenceOn_times to feedback_times
%
% CONTINUOUS REGRESSOR
% --------------------
%
%   wheelSpeed:
%       unsigned time derivative of wheel position, interpolated to
%       frameTimes
%
% OPTIONS
% -------
%
% lagWindow
%   Default lag range for discrete kernels, in seconds.
%   Default: [0 1.2].
%
% goCueLagWindow
%   Lag range for the goCue kernel relative to choiceMovement.
%   Default: [-0.4 0.6].
%
% choiceLagWindow
%   Lag range for the signed choice kernel relative to choiceMovement.
%   Default: [-1 0.2].
%
% actionLagWindow
%   Lag range for the action kernel relative to choiceMovement.
%   Default: [-0.2 1].
%
% frameRate
%   Sampling rate used to define the discrete lag grids. By default this is
%   inferred from the median interval in frameTimes.
%
% includeIntercept
%   Include a constant column. Default: true.
%
% contrastExponent
%   Exponent applied to left and right stimulus contrasts. Default: 1.
%
% centerBlock
%   Subtract 0.5 from probabilityLeft. Default: true.
%
% centerWheelSpeed
%   Subtract the mean from wheel speed. Default: false.
%
% scaleWheelSpeed
%   Divide wheel speed by its standard deviation. Default: true.
%
% wheelSpeedSmoothing
%   Smoothing window applied to wheel speed after interpolation, in
%   seconds. Set to zero to disable smoothing. Default: 0.
%
% LAG CONVENTION
% --------------
%
% A positive lag means that the design-matrix value occurs after the event.
% For example, the +0.4-second goCue column contains go-cue impulses shifted
% 0.4 seconds later in time.
%
% Negative choice and action lags represent activity preceding the
% choiceMovement event.

arguments
    trialsT table
    frameTimes (:,1) double
    wheelPosition (:,1) double
    wheelTimestamps (:,1) double

    options.lagWindow (1,2) double = [0 1.2]
    options.goCueLagWindow (1,2) double = [-0.6 0.6]
    options.choiceLagWindow (1,2) double = [-1 0.2]
    options.actionLagWindow (1,2) double = [-0.2 1]

    options.frameRate (1,1) double = NaN
    options.includeIntercept (1,1) logical = true
    options.contrastExponent (1,1) double {mustBePositive} = 1

    options.centerBlock (1,1) logical = true

    options.centerWheelSpeed (1,1) logical = false
    options.scaleWheelSpeed (1,1) logical = true
    options.wheelSpeedSmoothing (1,1) double = 0
end

%% Standardize input shapes

frameTimes = frameTimes(:);
wheelPosition = wheelPosition(:);
wheelTimestamps = wheelTimestamps(:);

%% Validate inputs

if isempty(frameTimes)
    error("buildIBLDesignMatrix:EmptyFrameTimes", ...
        "frameTimes cannot be empty.");
end

if numel(frameTimes) < 2
    error("buildIBLDesignMatrix:InsufficientFrameTimes", ...
        "At least two imaging frame times are required.");
end

if any(~isfinite(frameTimes))
    error("buildIBLDesignMatrix:InvalidFrameTimes", ...
        "frameTimes must contain only finite values.");
end

if any(diff(frameTimes) <= 0)
    error("buildIBLDesignMatrix:UnsortedFrameTimes", ...
        "frameTimes must be strictly increasing.");
end

if numel(wheelPosition) ~= numel(wheelTimestamps)
    error("buildIBLDesignMatrix:WheelSizeMismatch", ...
        "wheelPosition and wheelTimestamps must have equal lengths.");
end

validateLagWindow(options.goCueLagWindow, "goCueLagWindow");
validateLagWindow(options.lagWindow, "lagWindow");
validateLagWindow(options.choiceLagWindow, "choiceLagWindow");
validateLagWindow(options.actionLagWindow, "actionLagWindow");

if options.wheelSpeedSmoothing < 0
    error("buildIBLDesignMatrix:InvalidSmoothing", ...
        "wheelSpeedSmoothing must be nonnegative.");
end

requiredVariables = [
    "goCue_times"
    "contrastLeft"
    "contrastRight"
    "choiceMovement_times"
    "choice"
    "feedback_times"
    "feedbackType"
    "probabilityLeft"
    "quiescenceOn_times"
    ];

trialVariableNames = string(trialsT.Properties.VariableNames);
missingVariables = setdiff(requiredVariables, trialVariableNames);

if ~isempty(missingVariables)
    error("buildIBLDesignMatrix:MissingTrialVariables", ...
        "The trials table is missing: %s", ...
        strjoin(missingVariables, ", "));
end

T = numel(frameTimes);

%% Determine the kernel sampling rate

medianFrameInterval = median(diff(frameTimes), "omitnan");

if isnan(options.frameRate)
    frameRate = 1 / medianFrameInterval;
else
    if options.frameRate <= 0
        error("buildIBLDesignMatrix:InvalidFrameRate", ...
            "frameRate must be positive.");
    end

    frameRate = options.frameRate;
end

%% Define kernel-specific lag grids

[defaultLagFrames, defaultLagSeconds] = defineLagGrid( ...
    options.lagWindow, frameRate);

[goCueLagFrames, goCueLagSeconds] = defineLagGrid( ...
    options.goCueLagWindow, frameRate);

[choiceLagFrames, choiceLagSeconds] = defineLagGrid( ...
    options.choiceLagWindow, frameRate);

[actionLagFrames, actionLagSeconds] = defineLagGrid( ...
    options.actionLagWindow, frameRate);

fprintf("Imaging frame rate used for kernels: %.3f Hz\n", frameRate);
fprintf("Median observed frame interval: %.4f seconds\n", ...
    medianFrameInterval);

fprintf( ...
    "Default kernel range: %.3f to %.3f seconds (%d columns)\n", ...
    defaultLagSeconds(1), ...
    defaultLagSeconds(end), ...
    numel(defaultLagSeconds));

fprintf( ...
    "Go-cue kernel range: %.3f to %.3f seconds (%d columns)\n", ...
    goCueLagSeconds(1), ...
    goCueLagSeconds(end), ...
    numel(goCueLagSeconds));

fprintf( ...
    "Choice kernel range: %.3f to %.3f seconds (%d columns)\n", ...
    choiceLagSeconds(1), ...
    choiceLagSeconds(end), ...
    numel(choiceLagSeconds));

fprintf( ...
    "Action kernel range: %.3f to %.3f seconds (%d columns)\n", ...
    actionLagSeconds(1), ...
    actionLagSeconds(end), ...
    numel(actionLagSeconds));

%% Left- and right-stimulus event vectors

contrastLeft = double(trialsT.contrastLeft);
contrastRight = double(trialsT.contrastRight);

% IBL represents the absent stimulus on a given side as NaN.
contrastLeft(~isfinite(contrastLeft)) = 0;
contrastRight(~isfinite(contrastRight)) = 0;

% Nonlinear contrast transformation.
contrastLeft = contrastLeft .^ options.contrastExponent;
contrastRight = contrastRight .^ options.contrastExponent;

stimulusLeftImpulse = eventsToFrameVector( ...
    trialsT.goCue_times, ...
    contrastLeft, ...
    frameTimes);

stimulusRightImpulse = eventsToFrameVector( ...
    trialsT.goCue_times, ...
    contrastRight, ...
    frameTimes);

%% Go-cue event vector

goCueAmplitude = ones(height(trialsT), 1);

goCueImpulse = eventsToFrameVector( ...
    trialsT.goCue_times, ...
    goCueAmplitude, ...
    frameTimes);

%% Choice event vector

choiceAmplitude = double(trialsT.choice);

choiceImpulse = eventsToFrameVector( ...
    trialsT.choiceMovement_times, ...
    choiceAmplitude, ...
    frameTimes);

%% Action event vector

actionAmplitude = ones(height(trialsT), 1);

actionImpulse = eventsToFrameVector( ...
    trialsT.choiceMovement_times, ...
    actionAmplitude, ...
    frameTimes);

%% Feedback event vector

feedbackAmplitude = double(trialsT.feedbackType);

feedbackImpulse = eventsToFrameVector( ...
    trialsT.feedback_times, ...
    feedbackAmplitude, ...
    frameTimes);

%% Block-probability epoch regressor

blockAmplitude = double(trialsT.probabilityLeft);

if options.centerBlock
    % probabilityLeft values:
    %
    %   0.2 -> -0.3
    %   0.5 ->  0.0
    %   0.8 -> +0.3
    blockAmplitude = blockAmplitude - 0.5;
end

blockRegressor = epochsToFrameVector( ...
    trialsT.quiescenceOn_times, ...
    trialsT.feedback_times, ...
    blockAmplitude, ...
    frameTimes);

includeBlock = numel(unique(blockRegressor)) == 3;

%% Calculate wheel speed and align it to imaging frames

wheelSpeed = calculateWheelSpeed( ...
    wheelPosition, ...
    wheelTimestamps, ...
    frameTimes);

if options.wheelSpeedSmoothing > 0
    smoothingFrames = max( ...
        1, round(options.wheelSpeedSmoothing * frameRate));

    wheelSpeed = smoothdata( ...
        wheelSpeed, ...
        "movmean", ...
        smoothingFrames);
end

if options.centerWheelSpeed
    wheelSpeed = wheelSpeed - mean(wheelSpeed, "omitnan");
end

if options.scaleWheelSpeed
    wheelSpeedStandardDeviation = std(wheelSpeed, 0, "omitnan");

    if isfinite(wheelSpeedStandardDeviation) && ...
            wheelSpeedStandardDeviation > 0

        wheelSpeed = wheelSpeed ./ wheelSpeedStandardDeviation;
    end
end

%% Construct Toeplitz-style lag matrices

stimulusLeftLagged = makeLaggedEventMatrix( ...
    stimulusLeftImpulse, defaultLagFrames);

stimulusRightLagged = makeLaggedEventMatrix( ...
    stimulusRightImpulse, defaultLagFrames);

goCueLagged = makeLaggedEventMatrix( ...
    goCueImpulse, goCueLagFrames);

choiceLagged = makeLaggedEventMatrix( ...
    choiceImpulse, choiceLagFrames);

actionLagged = makeLaggedEventMatrix( ...
    actionImpulse, actionLagFrames);

feedbackLagged = makeLaggedEventMatrix( ...
    feedbackImpulse, defaultLagFrames);

%% Assemble the design matrix

X = zeros(T, 0);

predictorName = strings(0,1);
predictorType = strings(0,1);
predictorLag = zeros(0,1);
columnIndex = zeros(0,1);

if options.includeIntercept
    X(:,end + 1) = 1;

    predictorName(end + 1,1) = "intercept";
    predictorType(end + 1,1) = "constant";
    predictorLag(end + 1,1) = NaN;
    columnIndex(end + 1,1) = size(X,2);
end

[X, predictorName, predictorType, predictorLag, columnIndex] = ...
    appendLaggedPredictor( ...
        X, ...
        predictorName, ...
        predictorType, ...
        predictorLag, ...
        columnIndex, ...
        goCueLagged, ...
        "goCue", ...
        goCueLagSeconds);

[X, predictorName, predictorType, predictorLag, columnIndex] = ...
    appendLaggedPredictor( ...
        X, ...
        predictorName, ...
        predictorType, ...
        predictorLag, ...
        columnIndex, ...
        stimulusLeftLagged, ...
        "stimulusLeft", ...
        defaultLagSeconds);

[X, predictorName, predictorType, predictorLag, columnIndex] = ...
    appendLaggedPredictor( ...
        X, ...
        predictorName, ...
        predictorType, ...
        predictorLag, ...
        columnIndex, ...
        stimulusRightLagged, ...
        "stimulusRight", ...
        defaultLagSeconds);

[X, predictorName, predictorType, predictorLag, columnIndex] = ...
    appendLaggedPredictor( ...
        X, ...
        predictorName, ...
        predictorType, ...
        predictorLag, ...
        columnIndex, ...
        choiceLagged, ...
        "choice", ...
        choiceLagSeconds);

[X, predictorName, predictorType, predictorLag, columnIndex] = ...
    appendLaggedPredictor( ...
        X, ...
        predictorName, ...
        predictorType, ...
        predictorLag, ...
        columnIndex, ...
        actionLagged, ...
        "action", ...
        actionLagSeconds);

[X, predictorName, predictorType, predictorLag, columnIndex] = ...
    appendLaggedPredictor( ...
        X, ...
        predictorName, ...
        predictorType, ...
        predictorLag, ...
        columnIndex, ...
        feedbackLagged, ...
        "feedback", ...
        defaultLagSeconds);

if includeBlock
    [X, predictorName, predictorType, predictorLag, columnIndex] = ...
        appendSinglePredictor( ...
            X, ...
            predictorName, ...
            predictorType, ...
            predictorLag, ...
            columnIndex, ...
            blockRegressor, ...
            "block", ...
            "epoch", ...
            NaN);
end

[X, predictorName, predictorType, predictorLag, columnIndex] = ...
    appendSinglePredictor( ...
        X, ...
        predictorName, ...
        predictorType, ...
        predictorLag, ...
        columnIndex, ...
        wheelSpeed, ...
        "wheelSpeed", ...
        "continuous", ...
        0);

%% Predictor metadata

predictorInfo = table( ...
    columnIndex, ...
    predictorName, ...
    predictorType, ...
    predictorLag, ...
    'VariableNames', { ...
        'column', ...
        'name', ...
        'type', ...
        'lagSeconds'});

%% Return frame-aligned predictors for inspection

predictors = struct;

predictors.frameTimes = frameTimes;
predictors.frameRate = frameRate;
predictors.contrastExponent = options.contrastExponent;

% Frame-aligned event impulses.
predictors.stimulusLeftImpulse = stimulusLeftImpulse;
predictors.stimulusRightImpulse = stimulusRightImpulse;
predictors.goCueImpulse = goCueImpulse;
predictors.choiceImpulse = choiceImpulse;
predictors.actionImpulse = actionImpulse;
predictors.feedbackImpulse = feedbackImpulse;

% Non-kernel predictors.
if includeBlock
    predictors.block = blockRegressor;
end

predictors.wheelSpeed = wheelSpeed;

% Kernel-specific lag definitions.
predictors.defaultLagFrames = defaultLagFrames;
predictors.defaultLagSeconds = defaultLagSeconds;

predictors.goCueLagFrames = goCueLagFrames;
predictors.goCueLagSeconds = goCueLagSeconds;

predictors.choiceLagFrames = choiceLagFrames;
predictors.choiceLagSeconds = choiceLagSeconds;

predictors.actionLagFrames = actionLagFrames;
predictors.actionLagSeconds = actionLagSeconds;

% Retain generic lag fields for backwards compatibility.
predictors.lagFrames = defaultLagFrames;
predictors.lagSeconds = defaultLagSeconds;

% Expanded kernel matrices.
predictors.stimulusLeftLagged = stimulusLeftLagged;
predictors.stimulusRightLagged = stimulusRightLagged;
predictors.goCueLagged = goCueLagged;
predictors.choiceLagged = choiceLagged;
predictors.actionLagged = actionLagged;
predictors.feedbackLagged = feedbackLagged;

%% Report possible missing data

nTrials = height(trialsT);

nValidQuiescence = sum(isfinite(trialsT.quiescenceOn_times));
nValidStimulus = sum(isfinite(trialsT.goCue_times));
nValidMovement = sum(isfinite(trialsT.choiceMovement_times));
nValidFeedback = sum(isfinite(trialsT.feedback_times));

if nValidQuiescence < nTrials
    warning("buildIBLDesignMatrix:MissingQuiescenceTimes", ...
        "%d of %d trials have missing quiescenceOn_times.", ...
        nTrials - nValidQuiescence, nTrials);
end

if nValidStimulus < nTrials
    warning("buildIBLDesignMatrix:MissingStimulusTimes", ...
        "%d of %d trials have missing goCue_times.", ...
        nTrials - nValidStimulus, nTrials);
end

if nValidMovement < nTrials
    warning("buildIBLDesignMatrix:MissingMovementTimes", ...
        "%d of %d trials have missing choiceMovement_times.", ...
        nTrials - nValidMovement, nTrials);
end

if nValidFeedback < nTrials
    warning("buildIBLDesignMatrix:MissingFeedbackTimes", ...
        "%d of %d trials have missing feedback_times.", ...
        nTrials - nValidFeedback, nTrials);
end

end


function validateLagWindow(lagWindow, optionName)
%VALIDATELAGWINDOW Validate a two-element kernel lag window.

if any(~isfinite(lagWindow))
    error("buildIBLDesignMatrix:InvalidLagWindow", ...
        "%s must contain two finite values.", optionName);
end

if lagWindow(1) > lagWindow(2)
    error("buildIBLDesignMatrix:InvalidLagWindow", ...
        "The first value of %s must not exceed the second.", optionName);
end

end


function [lagFrames, lagSeconds] = defineLagGrid(lagWindow, frameRate)
%DEFINELAGGRID Convert a lag window in seconds into an integer frame grid.
%
% The requested limits are rounded to the nearest frame. The returned
% lagSeconds values therefore describe the exact lags represented by the
% design-matrix columns.

minimumLagFrames = round(lagWindow(1) * frameRate);
maximumLagFrames = round(lagWindow(2) * frameRate);

lagFrames = minimumLagFrames:maximumLagFrames;
lagSeconds = lagFrames ./ frameRate;

end

function eventVector = eventsToFrameVector( ...
    eventTimes, eventAmplitudes, frameTimes)
%EVENTSTOFRAMEVECTOR Map discrete events to their nearest imaging frames.
%
% Events outside the imaging interval are ignored. If multiple events are
% assigned to the same frame, their amplitudes are summed.

eventTimes = double(eventTimes(:));
eventAmplitudes = double(eventAmplitudes(:));
frameTimes = frameTimes(:);

if numel(eventTimes) ~= numel(eventAmplitudes)
    error("eventsToFrameVector:SizeMismatch", ...
        "eventTimes and eventAmplitudes must have equal lengths.");
end

T = numel(frameTimes);
eventVector = zeros(T,1);

validEvents = ...
    isfinite(eventTimes) & ...
    isfinite(eventAmplitudes) & ...
    eventTimes >= frameTimes(1) & ...
    eventTimes <= frameTimes(end);

eventTimes = eventTimes(validEvents);
eventAmplitudes = eventAmplitudes(validEvents);

if isempty(eventTimes)
    return
end

frameIndices = interp1( ...
    frameTimes, ...
    (1:T)', ...
    eventTimes, ...
    "nearest", ...
    NaN);

validIndices = isfinite(frameIndices);

frameIndices = round(frameIndices(validIndices));
eventAmplitudes = eventAmplitudes(validIndices);

eventVector = accumarray( ...
    frameIndices, ...
    eventAmplitudes, ...
    [T 1], ...
    @sum, ...
    0);

end


function epochVector = epochsToFrameVector( ...
    startTimes, endTimes, epochAmplitudes, frameTimes)
%EPOCHSTOFRAMEVECTOR Construct a frame-aligned epoch regressor.
%
% For each trial, the relevant imaging frames are filled with the trial's
% epoch amplitude. Invalid epochs are skipped.
%
% If epochs overlap, later rows in the trials table overwrite earlier rows.

startTimes = double(startTimes(:));
endTimes = double(endTimes(:));
epochAmplitudes = double(epochAmplitudes(:));
frameTimes = frameTimes(:);

if numel(startTimes) ~= numel(endTimes) || ...
        numel(startTimes) ~= numel(epochAmplitudes)

    error("epochsToFrameVector:SizeMismatch", ...
        "Start times, end times, and amplitudes must have equal lengths.");
end

epochVector = zeros(numel(frameTimes),1);

validEpochs = ...
    isfinite(startTimes) & ...
    isfinite(endTimes) & ...
    isfinite(epochAmplitudes) & ...
    endTimes >= startTimes;

validEpochIndices = find(validEpochs);

for index = 1:numel(validEpochIndices)
    trialIndex = validEpochIndices(index);

    firstFrame = find( ...
        frameTimes >= startTimes(trialIndex), 1, "first");

    lastFrame = find( ...
        frameTimes <= endTimes(trialIndex), 1, "last");

    if isempty(firstFrame) || isempty(lastFrame) || firstFrame > lastFrame
        continue
    end

    epochVector(firstFrame:lastFrame) = epochAmplitudes(trialIndex);
end

end


function wheelSpeed = calculateWheelSpeed( ...
    wheelPosition, wheelTimestamps, frameTimes)
%CALCULATEWHEELSPEED Compute unsigned wheel speed at imaging-frame times.
%
% Wheel position is cleaned and sorted, then interpolated onto the imaging
% frame time base. The signed derivative is calculated after interpolation,
% and its absolute value is returned:
%
%     wheelSpeed = abs(d(wheelPosition) / dt)
%
% Therefore, clockwise and counterclockwise wheel movements both produce
% positive speed values. Frames outside the recorded wheel interval are
% assigned zero.

wheelPosition = double(wheelPosition(:));
wheelTimestamps = double(wheelTimestamps(:));
frameTimes = double(frameTimes(:));

validSamples = ...
    isfinite(wheelPosition) & ...
    isfinite(wheelTimestamps);

wheelPosition = wheelPosition(validSamples);
wheelTimestamps = wheelTimestamps(validSamples);

if numel(wheelTimestamps) < 2
    error("calculateWheelSpeed:InsufficientSamples", ...
        "At least two valid wheel samples are required.");
end

% Sort samples by timestamp.
[wheelTimestamps, sortIndex] = sort(wheelTimestamps);
wheelPosition = wheelPosition(sortIndex);

% Average wheel positions at duplicate timestamps.
[uniqueTimestamps, ~, timestampGroup] = unique(wheelTimestamps);

if numel(uniqueTimestamps) < numel(wheelTimestamps)
    wheelPosition = accumarray( ...
        timestampGroup, ...
        wheelPosition, ...
        [], ...
        @mean);

    wheelTimestamps = uniqueTimestamps;
end

if numel(wheelTimestamps) < 2
    error("calculateWheelSpeed:InsufficientUniqueSamples", ...
        "At least two unique wheel timestamps are required.");
end

% Interpolate wheel position onto the imaging-frame time base.
wheelPositionAtFrames = interp1( ...
    wheelTimestamps, ...
    wheelPosition, ...
    frameTimes, ...
    "linear", ...
    NaN);

validFrameRange = ...
    frameTimes >= wheelTimestamps(1) & ...
    frameTimes <= wheelTimestamps(end);

wheelSpeed = zeros(size(frameTimes));

if sum(validFrameRange) >= 2
    wheelVelocity = gradient( ...
        wheelPositionAtFrames(validFrameRange), ...
        frameTimes(validFrameRange));

    wheelSpeed(validFrameRange) = abs(wheelVelocity);
end

% Frames outside the recorded wheel interval remain zero.
wheelSpeed(~isfinite(wheelSpeed)) = 0;

end

function laggedMatrix = makeLaggedEventMatrix( ...
    eventVector, lagFrames)
%MAKELAGGEDEVENTMATRIX Construct a Toeplitz-style event matrix.
%
% laggedMatrix(t,k) = eventVector(t - lagFrames(k)).
%
% A positive lag therefore shifts an event toward later imaging frames.
% A negative lag shifts an event toward earlier imaging frames.

eventVector = eventVector(:);
lagFrames = lagFrames(:)';

T = numel(eventVector);
nLags = numel(lagFrames);

laggedMatrix = zeros(T, nLags);

for lagIndex = 1:nLags
    shift = lagFrames(lagIndex);

    if shift == 0
        laggedMatrix(:,lagIndex) = eventVector;

    elseif shift > 0
        if shift < T
            laggedMatrix((1 + shift):T, lagIndex) = ...
                eventVector(1:(T - shift));
        end

    else
        advance = -shift;

        if advance < T
            laggedMatrix(1:(T - advance), lagIndex) = ...
                eventVector((1 + advance):T);
        end
    end
end

end


function [X, names, types, lags, columns] = appendLaggedPredictor( ...
    X, names, types, lags, columns, ...
    newPredictors, baseName, lagSeconds)
%APPENDLAGGEDPREDICTOR Append a complete event-kernel matrix.

if size(newPredictors,1) ~= size(X,1)
    error("appendLaggedPredictor:RowMismatch", ...
        "New predictors must have the same number of rows as X.");
end

if size(newPredictors,2) ~= numel(lagSeconds)
    error("appendLaggedPredictor:LagMismatch", ...
        "The number of predictor columns must equal the number of lags.");
end

firstColumn = size(X,2) + 1;

X = [X, newPredictors];

lastColumn = size(X,2);
nNewColumns = size(newPredictors,2);

names = [
    names
    repmat(string(baseName), nNewColumns, 1)
    ];

types = [
    types
    repmat("discrete-lagged", nNewColumns, 1)
    ];

lags = [
    lags
    lagSeconds(:)
    ];

columns = [
    columns
    (firstColumn:lastColumn)'
    ];

end


function [X, names, types, lags, columns] = appendSinglePredictor( ...
    X, names, types, lags, columns, ...
    newPredictor, predictorName, predictorType, predictorLag)
%APPENDSINGLEPREDICTOR Append one unlagged predictor to the design matrix.

newPredictor = newPredictor(:);

if numel(newPredictor) ~= size(X,1)
    error("appendSinglePredictor:RowMismatch", ...
        "The new predictor must have the same number of rows as X.");
end

X(:,end + 1) = newPredictor;

names(end + 1,1) = string(predictorName);
types(end + 1,1) = string(predictorType);
lags(end + 1,1) = predictorLag;
columns(end + 1,1) = size(X,2);

end