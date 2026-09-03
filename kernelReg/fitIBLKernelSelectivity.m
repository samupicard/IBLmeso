function results = fitIBLKernelSelectivity( ...
    F, frameTimes, trialsT, wheelPosition, wheelTimestamps, rankR, options)
%FITIBLKERNELSELECTIVITY Estimate unique predictor selectivity.
%
% results = fitIBLKernelSelectivity( ...
%     F, frameTimes, trialsT, wheelPosition, wheelTimestamps, rankR)
%
% For each predictor:
%
%   1. Fit all other predictors using reduced-rank Poisson regression.
%   2. Generate out-of-fold predictions of neural activity.
%   3. Subtract those predictions from the measured activity.
%   4. Retain the positive residual component.
%   5. Fit the held-out predictor to those residual responses.
%   6. Calculate cross-validated Poisson deviance explained.
%
% A neuron is considered selective when the unique cross-validated
% deviance explained exceeds options.selectivityThreshold.
%
% IMPORTANT
% ---------
% The additive residual F - FhatOther may contain negative values, whereas
% fitReducedRankGLMNet requires a nonnegative Poisson response. This
% implementation therefore uses:
%
%       FresPositive = max(F - FhatOther, 0)
%
% This measures activity above that predicted by the other regressors but
% does not model negative residuals.
%
% INPUTS
% ------
% F
%   T-by-N nonnegative neural activity matrix.
%
% frameTimes
%   T-by-1 imaging-frame timestamps.
%
% trialsT
%   IBL trials table.
%
% wheelPosition
%   Wheel-position samples.
%
% wheelTimestamps
%   Timestamps corresponding to wheelPosition.
%
% rankR
%   Reduced rank.
%
% OPTIONS
% -------
% contrastExponent
%   Exponent applied to stimulus contrast. Values below 1 boost
%   intermediate contrasts. Default: 1.
%
% lagWindow
%   Discrete event lag window in seconds. Default: [-1 2].
%
% wheelSpeedSmoothing
%   Wheel-speed smoothing window in seconds. Default: 0.1.
%
% nFolds
%   Number of outer cross-validation folds. Default: 5.
%
% foldId
%   Optional T-by-1 outer-fold assignment.
%
% blockedFolds
%   Use contiguous temporal folds when foldId is empty. Default: true.
%
% selectivityThreshold
%   Unique deviance-explained threshold. Default: 0.01.
%
% alpha
%   Elastic-net mixing parameter. Default: 0.01.
%
% lambdaChoice
%   "lambda_min" or "lambda_1se". Default: "lambda_min".
%
% innerFolds
%   Number of glmnet folds used within each training set. Default: 5.
%
% centerPredictors
%   Center predictors in fitReducedRankGLMNet. Default: true.
%
% standardizePredictors
%   Standardize temporal basis functions within glmnet. Default: true.
%
% useParallel
%   Use parallel neuron fitting where supported. Default: false.
%
% storeFoldModels
%   Store every cross-validation model. This may consume substantial
%   memory. Default: false.
%
% fitFinalModels
%   Fit full-data models for each predictor after cross-validation.
%   Default: true.
%
% OUTPUT
% ------
% results.predictorNames
%   P-by-1 string array of tested predictor names.
%
% results.deUnique
%   N-by-P cross-validated unique Poisson deviance explained.
%
% results.deUniqueFraction
%   N-by-P unique deviance divided by the full model's cross-validated
%   deviance explained.
%
% results.selective
%   N-by-P logical selectivity matrix.
%
% results.full.devianceExplained
%   Cross-validated deviance explained by the full model.
%
% results.predictor(p).otherModel
%   Full-data model excluding predictor p, when fitFinalModels is true.
%
% results.predictor(p).testModel
%   Full-data model fitting predictor p to positive residual activity.
%
% results.predictor(p).oofOtherPrediction
%   Out-of-fold prediction from all predictors other than p.
%
% results.predictor(p).oofResidualPrediction
%   Out-of-fold prediction of positive residual activity from predictor p.

arguments
    F double
    frameTimes (:,1) double
    trialsT table
    wheelPosition (:,1) double
    wheelTimestamps (:,1) double
    rankR (1,1) double {mustBeInteger,mustBePositive}

    options.contrastExponent (1,1) double {mustBePositive} = 1
    options.lagWindow (1,2) double = [0 1.2]
    options.wheelSpeedSmoothing (1,1) double {mustBeNonnegative} = 0.1

    options.inTrialMask (1,1) logical = false

    options.nFolds (1,1) double {mustBeInteger,mustBePositive} = 5
    options.foldId double = []
    options.blockedFolds (1,1) logical = true

    options.selectivityThreshold (1,1) double {mustBeNonnegative} = 0.01

    options.alpha (1,1) double {mustBeNonnegative,mustBeLessThanOrEqual(options.alpha,1)} = 0.01
    options.lambdaChoice (1,1) string = "lambda_min"
    options.innerFolds (1,1) double {mustBeInteger,mustBePositive} = 5

    options.centerPredictors (1,1) logical = true
    options.standardizePredictors (1,1) logical = true
    options.useParallel (1,1) logical = true

    options.storeFoldModels (1,1) logical = false
    options.fitFinalModels (1,1) logical = true
    options.verbose (1,1) logical = true
end

%% Validate neural data

nFrames = numel(frameTimes);

if size(F,1) ~= nFrames
    error( ...
        "fitIBLKernelSelectivity:FrameCountMismatch", ...
        "F has %d rows, but frameTimes contains %d samples.", ...
        size(F,1),nFrames);
end

if any(F(:) < 0)
    error( ...
        "fitIBLKernelSelectivity:NegativeActivity", ...
        "F must be nonnegative for Poisson regression.");
end

nNeurons = size(F,2);

%% Construct the complete design matrix

% fitReducedRankGLMNet estimates its own neuron-specific intercept, so the
% design matrix should not contain an additional constant column.

[Xfull,predictorInfo,predictors] = buildIBLDesignMatrix( ...
    trialsT, ...
    frameTimes, ...
    wheelPosition, ...
    wheelTimestamps, ...
    contrastExponent=options.contrastExponent, ...
    lagWindow=options.lagWindow, ...
    wheelSpeedSmoothing=options.wheelSpeedSmoothing, ...
    includeIntercept=false);

if size(Xfull,1) ~= nFrames
    error( ...
        "fitIBLKernelSelectivity:DesignMatrixSize", ...
        "The design matrix and neural activity have different row counts.");
end

predictorNamePerColumn = string(predictorInfo.name);

% Exclude any intercept-like column if one was nevertheless returned.
isIntercept = strcmpi(predictorNamePerColumn,"intercept");

predictorNames = unique( ...
    predictorNamePerColumn(~isIntercept), ...
    "stable");

nPredictors = numel(predictorNames);

%% Define common valid observations

validRows = ...
    isfinite(frameTimes) & ...
    all(isfinite(Xfull),2) & ...
    all(isfinite(F),2);

if options.inTrialMask
    inTrialMask = makeTrialFrameMask(trialsT,frameTimes);
    validRows = validRows & inTrialMask;
end


if nnz(validRows) < options.nFolds
    error( ...
        "fitIBLKernelSelectivity:TooFewObservations", ...
        "Too few valid observations for %d cross-validation folds.", ...
        options.nFolds);
end

%% Construct outer cross-validation folds

if isempty(options.foldId)

    foldId = nan(nFrames,1);

    validIndices = find(validRows);
    nValid = numel(validIndices);

    if options.blockedFolds
        % Contiguous temporal folds.
        foldAssignment = ceil( ...
            (1:nValid)' * options.nFolds / nValid);
    else
        % Approximately balanced randomly interleaved folds.
        randomOrder = randperm(nValid);
        foldAssignment = zeros(nValid,1);
        foldAssignment(randomOrder) = ...
            mod((0:nValid-1)',options.nFolds) + 1;
    end

    foldId(validIndices) = foldAssignment;

else

    suppliedFoldId = options.foldId(:);

    if numel(suppliedFoldId) == nFrames
        foldId = suppliedFoldId;

    elseif numel(suppliedFoldId) == nnz(validRows)
        foldId = nan(nFrames,1);
        foldId(validRows) = suppliedFoldId;

    else
        error( ...
            "fitIBLKernelSelectivity:FoldIdSize", ...
            "foldId must contain one value per frame or valid frame.");
    end
end

foldValues = unique(foldId(validRows));
foldValues = foldValues(isfinite(foldValues));

if numel(foldValues) < 2
    error( ...
        "fitIBLKernelSelectivity:InvalidFolds", ...
        "At least two cross-validation folds are required.");
end

%% Fit full model cross-validation benchmark

if options.verbose
    fprintf("Cross-validating full model...\n");
end

fullCV = crossValidateReducedRankModel( ...
    Xfull, ...
    F, ...
    validRows, ...
    foldId, ...
    foldValues, ...
    rankR, ...
    options);

%% Allocate predictor results

deUnique = nan(nNeurons,nPredictors);
deUniqueFraction = nan(nNeurons,nPredictors);
selective = false(nNeurons,nPredictors);

emptyPredictorResult = struct( ...
    "name","", ...
    "testColumns",[], ...
    "otherColumns",[], ...
    "devianceExplained",[], ...
    "uniqueFraction",[], ...
    "selective",[], ...
    "oofOtherPrediction",[], ...
    "oofPositiveResidual",[], ...
    "oofResidualPrediction",[], ...
    "otherModel",[], ...
    "testModel",[], ...
    "foldOtherModels",[], ...
    "foldTestModels",[]);

predictorResults = repmat(emptyPredictorResult,nPredictors,1);

%% Test each predictor

for predictorIndex = 1:nPredictors

    predictorName = predictorNames(predictorIndex);

    if options.verbose
        fprintf( ...
            "Testing predictor %d/%d: %s\n", ...
            predictorIndex,nPredictors,predictorName);
    end

    testColumns = predictorNamePerColumn == predictorName;
    otherColumns = ~testColumns & ~isIntercept;

    Xtest = Xfull(:,testColumns);
    Xother = Xfull(:,otherColumns);

    if isempty(Xtest)
        warning( ...
            "fitIBLKernelSelectivity:NoTestColumns", ...
            "No design-matrix columns found for predictor %s.", ...
            predictorName);
        continue
    end

    if isempty(Xother)
        warning( ...
            "fitIBLKernelSelectivity:NoOtherColumns", ...
            "No remaining predictors after removing %s.", ...
            predictorName);
        continue
    end

    predictorCV = crossValidateResidualModel( ...
        Xother, ...
        Xtest, ...
        F, ...
        validRows, ...
        foldId, ...
        foldValues, ...
        rankR, ...
        options);

    thisDE = predictorCV.devianceExplained;

    thisFraction = thisDE ./ fullCV.devianceExplained;

    % Avoid unstable ratios when the full model explains no positive
    % deviance.
    invalidFraction = ...
        ~isfinite(thisFraction) | ...
        fullCV.devianceExplained <= 0;

    thisFraction(invalidFraction) = NaN;

    thisSelective = ...
        isfinite(thisDE) & ...
        thisDE > options.selectivityThreshold;

    deUnique(:,predictorIndex) = thisDE;
    deUniqueFraction(:,predictorIndex) = thisFraction;
    selective(:,predictorIndex) = thisSelective;

    predictorResults(predictorIndex).name = predictorName;
    predictorResults(predictorIndex).testColumns = find(testColumns);
    predictorResults(predictorIndex).otherColumns = find(otherColumns);
    predictorResults(predictorIndex).devianceExplained = thisDE;
    predictorResults(predictorIndex).uniqueFraction = thisFraction;
    predictorResults(predictorIndex).selective = thisSelective;

    predictorResults(predictorIndex).oofOtherPrediction = ...
        predictorCV.otherPrediction;

    predictorResults(predictorIndex).oofPositiveResidual = ...
        predictorCV.positiveResidual;

    predictorResults(predictorIndex).oofResidualPrediction = ...
        predictorCV.residualPrediction;

    if options.storeFoldModels
        predictorResults(predictorIndex).foldOtherModels = ...
            predictorCV.otherModels;

        predictorResults(predictorIndex).foldTestModels = ...
            predictorCV.testModels;
    end

    %% Fit final full-data models for interpretation

    if options.fitFinalModels

        XotherValid = Xother(validRows,:);
        XtestValid = Xtest(validRows,:);
        Fvalid = F(validRows,:);

        otherModel = fitModel( ...
            XotherValid,Fvalid,rankR,options);

        FhatOther = predictModel(otherModel,XotherValid);

        positiveResidual = max(Fvalid - FhatOther,0);

        testModel = fitModel( ...
            XtestValid,positiveResidual,rankR,options);

        predictorResults(predictorIndex).otherModel = otherModel;
        predictorResults(predictorIndex).testModel = testModel;
    end
end

%% Fit final full model

if options.fitFinalModels

    if options.verbose
        fprintf("Fitting final full model...\n");
    end

    fullModel = fitModel( ...
        Xfull(validRows,:), ...
        F(validRows,:), ...
        rankR, ...
        options);

else
    fullModel = [];
end

%% Assemble output

results = struct;

results.predictorNames = predictorNames;
results.deUnique = deUnique;
results.deUniqueFraction = deUniqueFraction;
results.selective = selective;
results.selectivityThreshold = options.selectivityThreshold;

results.full = struct;
results.full.model = fullModel;
results.full.oofPrediction = fullCV.prediction;
results.full.devianceExplained = fullCV.devianceExplained;
results.full.modelDeviance = fullCV.modelDeviance;
results.full.nullDeviance = fullCV.nullDeviance;

results.predictor = predictorResults;

results.Xfull = Xfull;
results.predictorInfo = predictorInfo;
results.predictors = predictors;

results.validRows = validRows;
results.foldId = foldId;
results.rank = rankR;
results.options = options;

results.residualDefinition = ...
    "max(F - predictionFromOtherPredictors, 0)";

end


function cv = crossValidateReducedRankModel( ...
    X, F, validRows, foldId, foldValues, rankR, options)
%CROSSVALIDATEREDUCEDRANKMODEL Generate out-of-fold model predictions.

[nFrames,nNeurons] = size(F);

prediction = nan(nFrames,nNeurons);
nullPrediction = nan(nFrames,nNeurons);

if options.storeFoldModels
    models = cell(numel(foldValues),1);
else
    models = {};
end

for foldIndex = 1:numel(foldValues)

    foldValue = foldValues(foldIndex);

    trainRows = validRows & foldId ~= foldValue;
    testRows = validRows & foldId == foldValue;

    model = fitModel( ...
        X(trainRows,:), ...
        F(trainRows,:), ...
        rankR, ...
        options);

    prediction(testRows,:) = ...
        predictModel(model,X(testRows,:));

    trainingMean = mean(F(trainRows,:),1);
    trainingMean = max(trainingMean,eps);

    nullPrediction(testRows,:) = ...
        repmat(trainingMean,nnz(testRows),1);

    if options.storeFoldModels
        models{foldIndex} = model;
    end
end

[modelDeviance,nullDeviance,devianceExplained] = ...
    calculateCrossValidatedDeviance( ...
        F(validRows,:), ...
        prediction(validRows,:), ...
        nullPrediction(validRows,:));

cv = struct;
cv.prediction = prediction;
cv.nullPrediction = nullPrediction;
cv.modelDeviance = modelDeviance;
cv.nullDeviance = nullDeviance;
cv.devianceExplained = devianceExplained;
cv.models = models;

end


function cv = crossValidateResidualModel( ...
    Xother, Xtest, F, validRows, foldId, foldValues, rankR, options)
%CROSSVALIDATERESIDUALMODEL Cross-validate the two-stage residual model.

[nFrames,nNeurons] = size(F);

otherPrediction = nan(nFrames,nNeurons);
positiveResidual = nan(nFrames,nNeurons);
residualPrediction = nan(nFrames,nNeurons);
residualNullPrediction = nan(nFrames,nNeurons);

if options.storeFoldModels
    otherModels = cell(numel(foldValues),1);
    testModels = cell(numel(foldValues),1);
else
    otherModels = {};
    testModels = {};
end

for foldIndex = 1:numel(foldValues)

    foldValue = foldValues(foldIndex);

    trainRows = validRows & foldId ~= foldValue;
    testRows = validRows & foldId == foldValue;

    %% Stage 1: fit all predictors other than the tested predictor

    otherModel = fitModel( ...
        Xother(trainRows,:), ...
        F(trainRows,:), ...
        rankR, ...
        options);

    FhatOtherTrain = predictModel( ...
        otherModel,Xother(trainRows,:));

    FhatOtherTest = predictModel( ...
        otherModel,Xother(testRows,:));

    otherPrediction(testRows,:) = FhatOtherTest;

    %% Construct nonnegative residual targets

    residualTrain = max(F(trainRows,:) - FhatOtherTrain,0);
    residualTest = max(F(testRows,:) - FhatOtherTest,0);

    positiveResidual(testRows,:) = residualTest;

    %% Stage 2: fit only the tested predictor to residual activity

    if all(residualTrain(:) == 0)

        residualPrediction(testRows,:) = 0;
        residualNullPrediction(testRows,:) = 0;

        if options.storeFoldModels
            otherModels{foldIndex} = otherModel;
            testModels{foldIndex} = [];
        end

        continue
    end

    testModel = fitModel( ...
        Xtest(trainRows,:), ...
        residualTrain, ...
        rankR, ...
        options);

    residualPrediction(testRows,:) = ...
        predictModel(testModel,Xtest(testRows,:));

    residualTrainingMean = mean(residualTrain,1);
    residualTrainingMean = max(residualTrainingMean,eps);

    residualNullPrediction(testRows,:) = ...
        repmat(residualTrainingMean,nnz(testRows),1);

    if options.storeFoldModels
        otherModels{foldIndex} = otherModel;
        testModels{foldIndex} = testModel;
    end
end

[modelDeviance,nullDeviance,devianceExplained] = ...
    calculateCrossValidatedDeviance( ...
        positiveResidual(validRows,:), ...
        residualPrediction(validRows,:), ...
        residualNullPrediction(validRows,:));

cv = struct;
cv.otherPrediction = otherPrediction;
cv.positiveResidual = positiveResidual;
cv.residualPrediction = residualPrediction;
cv.residualNullPrediction = residualNullPrediction;
cv.modelDeviance = modelDeviance;
cv.nullDeviance = nullDeviance;
cv.devianceExplained = devianceExplained;
cv.otherModels = otherModels;
cv.testModels = testModels;

end


function model = fitModel(X,F,rankR,options)
%FITMODEL Call fitReducedRankGLMNet using common options.

maximumRank = min([rankR,size(X,2),size(F,2)]);

if maximumRank < 1
    error( ...
        "fitIBLKernelSelectivity:InvalidRank", ...
        "The effective reduced rank is less than one.");
end

model = fitReducedRankGLMNet( ...
    X, ...
    F, ...
    maximumRank, ...
    alpha=options.alpha, ...
    ...%relativeLambdaGrid = 10.^[-3 -2.25 -1.5 -0.75 0], ...
    relativeLambdaGrid = 10.^[-3 -1.5 0], ...
    validationFraction = 0.20, ...
    ridgeB = 1e-4, ...
    logOffset = 1e-3, ...
    ...%lambdaChoice=options.lambdaChoice, ...
    ...%nFolds=min(options.innerFolds,size(X,1)), ...
    standardizePredictors=options.standardizePredictors, ...
    centerPredictors=options.centerPredictors, ...
    useParallel=options.useParallel, ...
    storeCVFits=false);

end


function Fhat = predictModel(model,X)
%PREDICTMODEL Predict responses from uncentered predictor values.

linearPredictor = ...
    X * model.K + model.rawPredictorIntercept;

% Match the numerical protection used in fitReducedRankGLMNet.
linearPredictor = min(linearPredictor,700);

Fhat = exp(linearPredictor);

end


function [modelDeviance,nullDeviance,devianceExplained] = ...
    calculateCrossValidatedDeviance(y,muModel,muNull)
%CALCULATECROSSVALIDATEDDEVIANCE Poisson deviance for each neuron.

muModel = max(muModel,eps);
muNull = max(muNull,eps);

modelDeviance = sum(poissonDevianceTerms(y,muModel),1);
nullDeviance = sum(poissonDevianceTerms(y,muNull),1);

devianceExplained = ...
    1 - modelDeviance ./ nullDeviance;

invalid = ...
    ~isfinite(devianceExplained) | ...
    nullDeviance <= eps;

devianceExplained(invalid) = NaN;

devianceExplained = devianceExplained(:);
modelDeviance = modelDeviance(:);
nullDeviance = nullDeviance(:);

end


function terms = poissonDevianceTerms(y,mu)
%POISSONDEVIANCETERMS Observation-wise Poisson deviance terms.

terms = 2 * (mu - y);

positive = y > 0;

terms(positive) = 2 * ( ...
    y(positive) .* log(y(positive) ./ mu(positive)) ...
    - (y(positive) - mu(positive)));

end