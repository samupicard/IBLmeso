function model = fitReducedRankGLMNet(P, F, rankR, options)
%FITREDUCEDRANKGLMNET Poisson reduced-rank basis with neuron-wise glmnet.
%
% model = fitReducedRankGLMNet(P, F, rankR)
%
% Fits the positive-response model
%
%     Fhat = exp((P - predictorMean) * B * W + intercept)
%
% where:
%
%     P               T-by-E predictor matrix
%     B               E-by-R shared kernel basis
%     W               R-by-N neuron-specific weights
%     temporalBasis   (P - predictorMean)*B, size T-by-R
%
% The shared basis B is initialized using reduced-rank regression in
% log-response space:
%
%     log(F + logOffset) approximately equals Pc * B * W0
%
% Given B, each neuron's W(:,n) and intercept(n) are then fitted using
% Poisson elastic-net regression with a log link through glmnet.
%
% This is a two-stage approximation. It does not jointly optimize B and W
% under a Poisson likelihood.
%
% INPUTS
% ------
% P
%   T-by-E predictor matrix.
%
% F
%   T-by-N nonnegative deconvolved firing-rate matrix.
%
% rankR
%   Desired reduced rank R.
%
% OPTIONS
% -------
% alpha
%   Elastic-net mixing parameter used when fitting W:
%
%       alpha = 1      lasso
%       alpha near 0   approximately ridge
%       0 < alpha < 1  elastic net
%
%   Default: 0.01.
%
% lambdaChoice
%   "lambda_min" or "lambda_1se".
%   Default: "lambda_min".
%
% nFolds
%   Number of cross-validation folds. Default: 5.
%
% foldId
%   Optional fold assignment.
%
%   This may have one value per original observation or one value per
%   retained valid observation.
%
% standardizePredictors
%   Whether glmnet standardizes temporal basis functions internally.
%   Default: true.
%
% centerPredictors
%   Center P before constructing the shared basis. Default: true.
%
% ridgeB
%   Ridge penalty used to initialize B in log-response space.
%   Default: 1e-6.
%
% logOffset
%   Positive offset used only for initialization:
%
%       Ylog = log(F + logOffset)
%
%   This does not affect the subsequent Poisson glmnet fit.
%   Default: 1e-3.
%
% exposure
%   Optional T-by-1 exposure vector.
%
%   For equal-duration imaging frames, leave this empty. If frames have
%   unequal effective durations, the Poisson model can include
%   log(exposure) as an offset.
%
%   Default: [].
%
% useParallel
%   Fit neurons with parfor. Default: false.
%
% storeCVFits
%   Store all cvglmnet objects. These may require substantial memory for
%   large neuron counts. Default: true.
%
% maximumLinearPredictor
%   Upper numerical limit applied before exponentiation. This prevents
%   overflow when constructing Fhat manually. Default: 700.
%
% OUTPUTS
% -------
% model.B
%   E-by-R shared kernel basis.
%
% model.W
%   R-by-N neuron weights estimated by Poisson glmnet.
%
% model.K
%   E-by-N linear-predictor kernel matrix B*W.
%
% model.temporalBasis
%   T-by-R matrix (P - predictorMean)*B.
%
% model.intercept
%   Intercepts expressed for centered predictors.
%
% model.rawPredictorIntercept
%   Equivalent intercepts for predictions made directly from uncentered P:
%
%       eta = P*model.K + model.rawPredictorIntercept
%
% model.linearPredictor
%   T-by-N fitted log rates.
%
% model.Fhat
%   T-by-N positive fitted mean responses.
%
% model.deviance
%   Poisson deviance for each neuron.
%
% model.nullDeviance
%   Poisson deviance of an intercept-only model for each neuron.
%
% model.devianceExplained
%   1 - deviance/nullDeviance for each neuron.
%
% model.rSquared
%   Conventional squared-error R-squared on the response scale. This is
%   retained as a descriptive metric, but deviance explained is generally
%   more appropriate for Poisson regression.
%
% model.lambda
%   Selected lambda for each neuron.
%
% model.cvFit
%   Cell array of cvglmnet output objects if storeCVFits is true.

arguments
    P (:,:) double
    F (:,:) double
    rankR (1,1) double {mustBeInteger, mustBePositive}

    options.alpha (1,1) double = 0.01
    options.lambdaChoice (1,1) string = "lambda_min"
    options.nFolds (1,1) double {mustBeInteger, mustBePositive} = 5
    options.foldId (:,1) double = []
    % When empty, use the original cvglmnet procedure.
    %
    % When nonempty, tune each neuron over these fractions of its
    % neuron-specific lambdaMax using one contiguous validation block.
    options.relativeLambdaGrid (1,:) double = []

    % Fraction of retained time points used as the contiguous validation
    % block in relative-lambda-grid mode.
    options.validationFraction (1,1) double = 0.20

    options.standardizePredictors (1,1) logical = true
    options.centerPredictors (1,1) logical = true

    options.ridgeB (1,1) double {mustBeNonnegative} = 1e-6
    options.logOffset (1,1) double {mustBePositive} = 1e-3

    options.exposure (:,1) double = []

    options.useParallel (1,1) logical = false
    options.storeCVFits (1,1) logical = false

    options.maximumLinearPredictor (1,1) double = 700
end

%% Validate dimensions

[Toriginal, E] = size(P);
[responseT, N] = size(F);

if responseT ~= Toriginal
    error("fitReducedRankGLMNet:TimeMismatch", ...
        "P and F must have the same number of rows.");
end

if rankR > min([E, N, Toriginal])
    error("fitReducedRankGLMNet:RankTooLarge", ...
        "rankR cannot exceed min(E,N,T).");
end

if options.alpha < 0 || options.alpha > 1
    error("fitReducedRankGLMNet:InvalidAlpha", ...
        "alpha must lie between zero and one.");
end

if ~ismember(options.lambdaChoice, ["lambda_min", "lambda_1se"])
    error("fitReducedRankGLMNet:InvalidLambdaChoice", ...
        "lambdaChoice must be lambda_min or lambda_1se.");
end

if options.maximumLinearPredictor <= 0
    error("fitReducedRankGLMNet:InvalidMaximumLinearPredictor", ...
        "maximumLinearPredictor must be positive.");
end

if ~isempty(options.relativeLambdaGrid)
    if any(~isfinite(options.relativeLambdaGrid)) || ...
            any(options.relativeLambdaGrid <= 0) || ...
            any(options.relativeLambdaGrid > 1)

        error("fitReducedRankGLMNet:InvalidRelativeLambdaGrid", ...
            ["relativeLambdaGrid must contain finite values in the " ...
             "interval (0,1]."]);
    end
end

if options.validationFraction <= 0 || ...
        options.validationFraction >= 1

    error("fitReducedRankGLMNet:InvalidValidationFraction", ...
        "validationFraction must lie strictly between zero and one.");
end

%% Validate nonnegative responses

finiteResponses = F(isfinite(F));

if any(finiteResponses < 0)
    minimumResponse = min(finiteResponses);

    error("fitReducedRankGLMNet:NegativeResponse", ...
        ["Poisson regression requires nonnegative responses. " ...
         "The minimum observed response is %.6g."], ...
        minimumResponse);
end

%% Validate exposure

if isempty(options.exposure)
    exposureOriginal = ones(Toriginal,1);
else
    if numel(options.exposure) ~= Toriginal
        error("fitReducedRankGLMNet:ExposureSizeMismatch", ...
            "exposure must contain one value per row of P.");
    end

    exposureOriginal = options.exposure(:);

    if any(~isfinite(exposureOriginal)) || any(exposureOriginal <= 0)
        error("fitReducedRankGLMNet:InvalidExposure", ...
            "All exposure values must be finite and strictly positive.");
    end
end

%% Remove rows containing invalid predictor or response values

validRows = ...
    all(isfinite(P), 2) & ...
    all(isfinite(F), 2) & ...
    isfinite(exposureOriginal) & ...
    exposureOriginal > 0;

if ~all(validRows)
    warning("fitReducedRankGLMNet:InvalidRows", ...
        "Removing %d rows containing invalid values.", ...
        sum(~validRows));
end

Pfit = P(validRows,:);
Ffit = F(validRows,:);
exposure = exposureOriginal(validRows);

T = size(Pfit,1);

if rankR > min([E, N, T])
    error("fitReducedRankGLMNet:RankTooLargeAfterFiltering", ...
        "After filtering invalid rows, rankR exceeds min(E,N,T).");
end

%% Center predictors

if options.centerPredictors
    predictorMean = mean(Pfit, 1);
else
    predictorMean = zeros(1,E);
end

Pc = Pfit - predictorMean;

%% Initialize the shared basis in log-response space

% We cannot center the raw response for Poisson regression.
%
% Instead, estimate the initial shared task subspace from:
%
%     Ylog = log(F + logOffset)
%
% The neuron-wise mean in log space is removed only for constructing B.
% The final Poisson model fits its own intercept for every neuron.

Ylog = log(Ffit + options.logOffset);

logResponseMean = mean(Ylog, 1);
YlogCentered = Ylog - logResponseMean;

% Ridge-stabilized multivariate regression:
%
%   KlogFull = argmin_K ||YlogCentered - Pc*K||_F^2
%                        + ridgeB*||K||_F^2

ridgeMatrix = options.ridgeB * eye(E);

KlogFull = (Pc' * Pc + ridgeMatrix) \ ...
    (Pc' * YlogCentered);

%% Obtain rank-R response-space directions

YlogPredictedFull = Pc * KlogFull;

[~, singularValues, V] = svd(YlogPredictedFull, "econ");

Vr = V(:,1:rankR);

%% Construct B and the temporal basis

% Initial reduced-rank log-kernel factorization:
%
%     B         = KlogFull * Vr
%     initialW  = Vr'
%
% The final W is refitted under a Poisson objective.

B = KlogFull * Vr;

temporalBasis = Pc * B;

%% Construct tuning partitions

useRelativeLambdaGrid = ~isempty(options.relativeLambdaGrid);

if useRelativeLambdaGrid
    % Use the final contiguous portion of the recording for validation.
    nValidation = max(1, round(options.validationFraction * T));

    if nValidation >= T
        error("fitReducedRankGLMNet:ValidationSetTooLarge", ...
            "The validation set leaves no observations for training.");
    end

    validationRows = false(T,1);
    validationRows((T - nValidation + 1):T) = true;
    trainingRows = ~validationRows;

    relativeLambdaGrid = unique( ...
        sort(options.relativeLambdaGrid(:)', "descend"), ...
        "stable");

    foldId = [];
    nFolds = NaN;

else
    relativeLambdaGrid = [];
    trainingRows = true(T,1);
    validationRows = false(T,1);

    if isempty(options.foldId)
        foldId = makeContiguousFolds(T, options.nFolds);
    else
        suppliedFoldId = options.foldId(:);

        if numel(suppliedFoldId) == Toriginal
            foldId = suppliedFoldId(validRows);

        elseif numel(suppliedFoldId) == T
            foldId = suppliedFoldId;

        else
            error("fitReducedRankGLMNet:FoldSizeMismatch", ...
                ["foldId must contain one value per original time point " ...
                 "or one value per retained valid time point."]);
        end
    end

    if any(~isfinite(foldId)) || any(foldId < 1)
        error("fitReducedRankGLMNet:InvalidFoldId", ...
            "foldId must contain finite positive fold labels.");
    end

    [~, ~, foldId] = unique(foldId, "stable");

    nFolds = numel(unique(foldId));

    if nFolds < 3
        error("fitReducedRankGLMNet:TooFewFolds", ...
            "glmnet cross-validation requires at least three folds.");
    end
end
%% Configure glmnet

glmOptions = glmnetSet;

glmOptions.alpha = options.alpha;
glmOptions.standardize = options.standardizePredictors;
glmOptions.intr = true;

% glmnet's Poisson implementation supports an offset through the options
% structure in some MATLAB wrapper versions. Because wrapper versions vary,
% exposure is incorporated into the response and prediction explicitly
% below only when exposure is constant.
%
% For ordinary fixed-duration imaging frames, exposure is all ones.

if any(abs(exposure - exposure(1)) > ...
        eps(max(abs(exposure))) * 10)

    warning("fitReducedRankGLMNet:VariableExposure", ...
        ["Variable exposure was supplied, but this function does not " ...
         "pass an offset into glmnet because MATLAB wrapper versions " ...
         "differ in offset support. Fit rates F./exposure or extend the " ...
         "glmnet call to include log(exposure) explicitly."]);
end

%% Fit neuron-specific Poisson models

W = zeros(rankR,N);
intercept = zeros(1,N);
selectedLambda = nan(1,N);
selectedRelativeLambda = nan(1,N);
lambdaMax = nan(1,N);

if useRelativeLambdaGrid
    validationDeviance = nan(N, numel(relativeLambdaGrid));
    cvFit = {};
else
    validationDeviance = [];
    selectedRelativeLambda = nan(1,N);
    lambdaMax = nan(1,N);

    if options.storeCVFits
        cvFit = cell(1,N);
    else
        cvFit = {};
    end
end

if useRelativeLambdaGrid
    if options.storeCVFits
        warning("fitReducedRankGLMNet:CVFitsUnavailable", ...
            ["storeCVFits is ignored in relative-lambda-grid mode " ...
             "because cvglmnet is not used."]);
    end

    if options.useParallel
        parfor neuronIndex = 1:N
            [W(:,neuronIndex), ...
                intercept(neuronIndex), ...
                selectedLambda(neuronIndex), ...
                selectedRelativeLambda(neuronIndex), ...
                lambdaMax(neuronIndex), ...
                validationDeviance(neuronIndex,:)] = ...
                fitSingleNeuronRelativeLambdaGrid( ...
                    temporalBasis, ...
                    Ffit(:,neuronIndex), ...
                    glmOptions, ...
                    trainingRows, ...
                    validationRows, ...
                    relativeLambdaGrid, ...
                    options.maximumLinearPredictor);
        end
    else
        for neuronIndex = 1:N
            [W(:,neuronIndex), ...
                intercept(neuronIndex), ...
                selectedLambda(neuronIndex), ...
                selectedRelativeLambda(neuronIndex), ...
                lambdaMax(neuronIndex), ...
                validationDeviance(neuronIndex,:)] = ...
                fitSingleNeuronRelativeLambdaGrid( ...
                    temporalBasis, ...
                    Ffit(:,neuronIndex), ...
                    glmOptions, ...
                    trainingRows, ...
                    validationRows, ...
                    relativeLambdaGrid, ...
                    options.maximumLinearPredictor);

            if mod(neuronIndex,100) == 0 || neuronIndex == N
                fprintf("Fitted %d of %d neurons.\n", neuronIndex, N);
            end
        end
    end

else
    if options.useParallel
        if options.storeCVFits
            parfor neuronIndex = 1:N
                [W(:,neuronIndex), ...
                    intercept(neuronIndex), ...
                    selectedLambda(neuronIndex), ...
                    cvFit{neuronIndex}] = fitSingleNeuronPoisson( ...
                        temporalBasis, ...
                        Ffit(:,neuronIndex), ...
                        glmOptions, ...
                        nFolds, ...
                        foldId, ...
                        options.lambdaChoice);
            end
        else
            parfor neuronIndex = 1:N
                [W(:,neuronIndex), ...
                    intercept(neuronIndex), ...
                    selectedLambda(neuronIndex)] = ...
                    fitSingleNeuronPoisson( ...
                        temporalBasis, ...
                        Ffit(:,neuronIndex), ...
                        glmOptions, ...
                        nFolds, ...
                        foldId, ...
                        options.lambdaChoice);
            end
        end
    else
        for neuronIndex = 1:N
            if options.storeCVFits
                [W(:,neuronIndex), ...
                    intercept(neuronIndex), ...
                    selectedLambda(neuronIndex), ...
                    cvFit{neuronIndex}] = fitSingleNeuronPoisson( ...
                        temporalBasis, ...
                        Ffit(:,neuronIndex), ...
                        glmOptions, ...
                        nFolds, ...
                        foldId, ...
                        options.lambdaChoice);
            else
                [W(:,neuronIndex), ...
                    intercept(neuronIndex), ...
                    selectedLambda(neuronIndex)] = ...
                    fitSingleNeuronPoisson( ...
                        temporalBasis, ...
                        Ffit(:,neuronIndex), ...
                        glmOptions, ...
                        nFolds, ...
                        foldId, ...
                        options.lambdaChoice);
            end

            if mod(neuronIndex,100) == 0 || neuronIndex == N
                fprintf("Fitted %d of %d neurons.\n", neuronIndex, N);
            end
        end
    end
end

%% Reconstruct kernels and fitted mean responses

% K operates in the linear-predictor/log-rate space.
K = B * W;

linearPredictorFit = temporalBasis * W + intercept;

% Prevent both numerical overflow and underflow during exponentiation.
% poissonDeviance requires fitted means to be finite and strictly positive.
minimumLinearPredictor = log(realmin("double"));

linearPredictorForExp = min( ...
    max(linearPredictorFit, minimumLinearPredictor), ...
    options.maximumLinearPredictor);

FhatFit = exp(linearPredictorForExp);

% Equivalent intercept when using the original, uncentered P:
%
%   (P - predictorMean)*K + intercept
%       = P*K + intercept - predictorMean*K

rawPredictorIntercept = intercept - predictorMean * K;

%% Return predictions aligned to original rows

linearPredictor = nan(Toriginal,N);
linearPredictor(validRows,:) = linearPredictorFit;

Fhat = nan(Toriginal,N);
Fhat(validRows,:) = FhatFit;

%% Calculate Poisson deviance

% Calculate deviance one neuron at a time so that neurons with an
% identically zero response can be handled explicitly. Their saturated
% model and intercept-only null model both have deviance zero, while
% deviance explained is undefined because it would require 0/0.
deviance = nan(1,N);
nullDeviance = nan(1,N);
nullMean = mean(Ffit, 1);
allZeroNeurons = all(Ffit == 0, 1);

for neuronIndex = 1:N
    response = Ffit(:,neuronIndex);

    if allZeroNeurons(neuronIndex)
        deviance(neuronIndex) = 0;
        nullDeviance(neuronIndex) = 0;
        continue
    end

    fittedMean = max( ...
        FhatFit(:,neuronIndex), ...
        realmin("double"));

    nullRate = max( ...
        nullMean(neuronIndex), ...
        realmin("double"));

    nullPrediction = repmat(nullRate, T, 1);

    deviance(neuronIndex) = poissonDeviance( ...
        response, ...
        fittedMean);

    nullDeviance(neuronIndex) = poissonDeviance( ...
        response, ...
        nullPrediction);
end

devianceExplained = nan(1,N);
validDevianceExplained = ...
    isfinite(deviance) & ...
    isfinite(nullDeviance) & ...
    nullDeviance > 0;

devianceExplained(validDevianceExplained) = ...
    1 - deviance(validDevianceExplained) ./ ...
        nullDeviance(validDevianceExplained);

populationDeviance = sum(deviance, "omitnan");
populationNullDeviance = sum(nullDeviance, "omitnan");

if populationNullDeviance > 0
    populationDevianceExplained = ...
        1 - populationDeviance / populationNullDeviance;
else
    populationDevianceExplained = NaN;
end

%% Also calculate response-scale squared-error metrics

residual = Ffit - FhatFit;

sumSquaredError = sum(residual.^2, 1);
totalSumSquares = sum((Ffit - mean(Ffit,1)).^2, 1);

rSquared = 1 - sumSquaredError ./ totalSumSquares;
rSquared(totalSumSquares == 0) = NaN;

populationSquaredError = sum(residual(:).^2);
populationTotalSumSquares = sum( ...
    (Ffit - mean(Ffit,1)).^2, ...
    "all");

if populationTotalSumSquares > 0
    populationR2 = 1 - ...
        populationSquaredError / populationTotalSumSquares;
else
    populationR2 = NaN;
end

%% Calculate mean Poisson log likelihood, ignoring constants

% The omitted term log(y!) does not depend on model parameters.
poissonLogLikelihood = sum( ...
    Ffit .* log(max(FhatFit, realmin)) - FhatFit, ...
    1);

%% Store outputs

model = struct;

model.family = "poisson";
model.link = "log";

model.B = B;
model.W = W;
model.K = K;

model.temporalBasis = temporalBasis;

model.intercept = intercept;
model.rawPredictorIntercept = rawPredictorIntercept;

model.linearPredictor = linearPredictor;
model.Fhat = Fhat;

model.predictorMean = predictorMean;
model.logResponseMean = logResponseMean;

model.rank = rankR;
model.alpha = options.alpha;
model.lambdaChoice = options.lambdaChoice;
model.lambda = selectedLambda;

if useRelativeLambdaGrid
    model.tuningMode = "relativeLambdaGrid";
else
    model.tuningMode = "crossValidation";
end

model.relativeLambdaGrid = relativeLambdaGrid;
model.selectedRelativeLambda = selectedRelativeLambda;
model.lambdaMax = lambdaMax;
model.validationDeviance = validationDeviance;
model.trainingRows = trainingRows;
model.validationRows = validationRows;
model.validationFraction = options.validationFraction;

model.foldId = foldId;
model.validRows = validRows;
model.allZeroNeurons = allZeroNeurons;
model.nullMean = nullMean;

model.deviance = deviance;
model.nullDeviance = nullDeviance;
model.devianceExplained = devianceExplained;
model.populationDevianceExplained = ...
    populationDevianceExplained;

model.rSquared = rSquared;
model.populationR2 = populationR2;

model.poissonLogLikelihood = poissonLogLikelihood;

model.exposure = exposureOriginal;

model.cvFit = cvFit;

% Initialization quantities.
model.fullLogKernel = KlogFull;
model.responseSubspace = Vr;
model.initialW = Vr';
model.initialReducedRankLogKernel = B * Vr';
model.singularValues = diag(singularValues);
model.logOffset = options.logOffset;

end


function [weights, intercept, lambda, cvFit] = ...
    fitSingleNeuronPoisson( ...
        temporalBasis, response, glmOptions, ...
        nFolds, foldId, lambdaChoice)
%FITSINGLENEURONPOISSON Fit one Poisson elastic-net model.
%
% glmnet models:
%
%     log(E[y|X]) = intercept + X*weights

if any(response < 0) || any(~isfinite(response))
    error("fitSingleNeuronPoisson:InvalidResponse", ...
        "Poisson responses must be finite and nonnegative.");
end

% A neuron that is identically zero cannot be fitted meaningfully with a
% finite Poisson intercept. Return a near-zero constant prediction.
if all(response == 0)
    weights = zeros(size(temporalBasis,2),1);
    intercept = log(realmin);
    lambda = NaN;
    cvFit = [];
    return
end

cvFit = cvglmnet( ...
    temporalBasis, ...
    response, ...
    "poisson", ...
    glmOptions, ...
    "deviance", ...
    nFolds, ...
    foldId, ...
    false);

coefficients = cvglmnetCoef( ...
    cvFit, ...
    char(lambdaChoice));

% The first coefficient is the intercept.
intercept = full(coefficients(1));
weights = full(coefficients(2:end));

switch lambdaChoice
    case "lambda_min"
        lambda = cvFit.lambda_min;

    case "lambda_1se"
        lambda = cvFit.lambda_1se;
end

end



function [weights, intercept, selectedLambda, ...
    selectedRelativeLambda, lambdaMaxFull, ...
    validationDeviance] = ...
    fitSingleNeuronRelativeLambdaGrid( ...
        temporalBasis, response, glmOptions, ...
        trainingRows, validationRows, ...
        relativeLambdaGrid, maximumLinearPredictor)
%FITSINGLENEURONRELATIVELAMBDAGRID
% Select a neuron-specific relative lambda using one contiguous validation
% block, then refit on all observations.
%
% Two glmnet paths are fitted:
%
%   1. A training-only path for selecting lambda/lambdaMax.
%   2. A full-data path from which the final coefficients are extracted.
%
% This is substantially cheaper than K-fold cvglmnet.

if any(response < 0) || any(~isfinite(response))
    error("fitSingleNeuronRelativeLambdaGrid:InvalidResponse", ...
        "Poisson responses must be finite and nonnegative.");
end

nPredictors = size(temporalBasis,2);
nCandidates = numel(relativeLambdaGrid);

validationDeviance = nan(1,nCandidates);

if all(response == 0)
    weights = zeros(nPredictors,1);
    intercept = log(realmin("double"));
    selectedLambda = NaN;
    selectedRelativeLambda = NaN;
    lambdaMaxFull = NaN;
    return
end

Xtrain = temporalBasis(trainingRows,:);
ytrain = response(trainingRows);

Xvalidation = temporalBasis(validationRows,:);
yvalidation = response(validationRows);

% There is no estimable training model when every training response is
% zero. Return a finite intercept-only model based on the complete response.
if all(ytrain == 0)
    weights = zeros(nPredictors,1);
    intercept = log(max(mean(response), realmin("double")));
    selectedLambda = NaN;
    selectedRelativeLambda = 1;
    lambdaMaxFull = NaN;
    return
end

%% Fit one training lambda path

trainingFit = glmnet( ...
    Xtrain, ...
    ytrain, ...
    "poisson", ...
    glmOptions);

trainingLambda = trainingFit.lambda(:);
lambdaMaxTraining = max(trainingLambda);

targetTrainingLambda = ...
    lambdaMaxTraining .* relativeLambdaGrid(:);

trainingIndices = nearestLambdaIndices( ...
    trainingLambda, ...
    targetTrainingLambda);

minimumLinearPredictor = log(realmin("double"));

%% Evaluate candidate relative lambdas

for candidateIndex = 1:nCandidates
    pathIndex = trainingIndices(candidateIndex);

    candidateIntercept = full(trainingFit.a0(pathIndex));
    candidateWeights = full(trainingFit.beta(:,pathIndex));

    etaValidation = ...
        Xvalidation * candidateWeights + candidateIntercept;

    etaValidation = min( ...
        max(etaValidation, minimumLinearPredictor), ...
        maximumLinearPredictor);

    predictionValidation = max( ...
        exp(etaValidation), ...
        realmin("double"));

    validationDeviance(candidateIndex) = ...
        poissonDeviance( ...
            yvalidation, ...
            predictionValidation);
end

finiteCandidates = isfinite(validationDeviance);

if any(finiteCandidates)
    finiteIndices = find(finiteCandidates);
    [~, localBestIndex] = min( ...
        validationDeviance(finiteCandidates));

    bestCandidateIndex = finiteIndices(localBestIndex);
else
    % Conservative fallback.
    bestCandidateIndex = find( ...
        relativeLambdaGrid == max(relativeLambdaGrid), ...
        1, ...
        "first");
end

selectedRelativeLambda = ...
    relativeLambdaGrid(bestCandidateIndex);

%% Refit one full-data lambda path

fullFit = glmnet( ...
    temporalBasis, ...
    response, ...
    "poisson", ...
    glmOptions);

fullLambda = fullFit.lambda(:);
lambdaMaxFull = max(fullLambda);

targetFullLambda = ...
    lambdaMaxFull * selectedRelativeLambda;

fullPathIndex = nearestLambdaIndices( ...
    fullLambda, ...
    targetFullLambda);

fullPathIndex = fullPathIndex(1);

selectedLambda = fullLambda(fullPathIndex);
intercept = full(fullFit.a0(fullPathIndex));
weights = full(fullFit.beta(:,fullPathIndex));

end


function indices = nearestLambdaIndices(lambdaPath, targetLambda)
%NEARESTLAMBDAINDICES Locate the closest path value in log-lambda space.

lambdaPath = lambdaPath(:);
targetLambda = targetLambda(:);

if any(~isfinite(lambdaPath)) || any(lambdaPath <= 0)
    error("nearestLambdaIndices:InvalidLambdaPath", ...
        "The glmnet lambda path must be finite and strictly positive.");
end

indices = zeros(numel(targetLambda),1);

logLambdaPath = log(lambdaPath);

for targetIndex = 1:numel(targetLambda)
    target = max(targetLambda(targetIndex), realmin("double"));

    [~, indices(targetIndex)] = min( ...
        abs(logLambdaPath - log(target)));
end

end
