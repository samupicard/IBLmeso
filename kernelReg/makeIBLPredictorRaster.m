function [predictorRaster, predictorLabels] = makeIBLPredictorRaster( ...
    predictors, options)
%MAKEIBLPREDICTORRASTER Create a compact signed task-predictor raster.
%
% [predictorRaster, predictorLabels] = makeIBLPredictorRaster(predictors)
%
% Output
%   predictorRaster
%       nPredictors-by-nFrames display matrix.
%
%   predictorLabels
%       Labels corresponding to the raster rows.
%
% Color/sign convention
%   Positive values -> red
%   Zero            -> white
%   Negative values -> blue
%
% Stimulus:
%   left contrast  -> positive
%   right contrast -> negative
%
% Choice:
%   left choice  -> positive
%   right choice -> negative
%
% This assumes predictors.choice uses the IBL convention:
%   -1 = left
%   +1 = right

arguments
    predictors struct
    options.includeBlock (1,1) logical = true
    options.includeWheelSpeed (1,1) logical = true
    options.normalizeRows (1,1) logical = true
end

%% Required discrete predictors

stimulusLeft = predictors.stimulusLeftImpulse(:).';
stimulusRight = predictors.stimulusRightImpulse(:).';
choice = predictors.choiceImpulse(:).';
feedback = predictors.feedbackImpulse(:).';

nFrames = numel(stimulusLeft);

assert(numel(stimulusRight) == nFrames, ...
    "stimulusLeft and stimulusRight have different lengths.");

assert(numel(choice) == nFrames, ...
    "choice has a different length from the stimulus predictors.");

assert(numel(feedback) == nFrames, ...
    "feedback has a different length from the stimulus predictors.");

%% Combine left and right stimulus into one signed predictor

% Positive/red = left stimulus
% Negative/blue = right stimulus
signedStimulus = stimulusLeft - stimulusRight;

%% Reorient choice so left is positive/red

% IBL coding:
%   choice = -1 means leftwards (CW)
%   choice = +1 means rightwards (CCW)

signedChoice = choice;

%% Assemble raster

predictorRaster = [
    signedStimulus
    signedChoice
    feedback
    ];

predictorLabels = [
    "Stimulus"
    "Choice"
    "Feedback"
    ];

%% Add block regressor

if options.includeBlock
    block = predictors.block(:).';

    assert(numel(block) == nFrames, ...
        "block has a different number of frames.");

    predictorRaster(end+1,:) = block;
    predictorLabels(end+1) = "Block";
end

%% Add wheel-speed regressor

if options.includeWheelSpeed
    wheelSpeed = predictors.wheelSpeed(:).';

    assert(numel(wheelSpeed) == nFrames, ...
        "wheelSpeed has a different number of frames.");

    predictorRaster(end+1,:) = wheelSpeed;
    predictorLabels(end+1) = "Wheel speed";
end

%% Replace nonfinite display values

predictorRaster(~isfinite(predictorRaster)) = 0;

%% Normalize each row independently for visualization

if options.normalizeRows
    maximumAbsoluteValue = max(abs(predictorRaster),[],2);

    nonzeroRows = maximumAbsoluteValue > 0;

    predictorRaster(nonzeroRows,:) = ...
        predictorRaster(nonzeroRows,:) ./ ...
        maximumAbsoluteValue(nonzeroRows);
end

end