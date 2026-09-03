function model = fitDirectGLMNet(P, F, options)
%FITDIRECTGLMNET Fit each neuron directly using Poisson elastic-net GLMs.
%
% model = fitDirectGLMNet(P, F)
%
% Fits independently for each neuron:
%
%     log(E[F(:,n) | P]) = P * K(:,n) + intercept(n)
%
% Therefore:
%
%     Fhat(:,n) = exp(P * K(:,n) + intercept(n))
%
% INPUTS
% ------
% P
%   T-by-E predictor matrix.
%
% F
%   T-by-N nonnegative neural response matrix.
%
% OPTIONS
% -------
% alpha
%   Elastic-net mixing parameter:
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
%   Optional T-by-1 fold assignment. If empty, contiguous temporal folds
%   are generated.
%
% standardizePredictors
%   Let glmnet standardize predictors internally. Default: true.
%
% includeGLMNetIntercept
%   Estimate an intercept in glmnet. Default: true.
%
% interceptColumn
%   Index of an explicit intercept column already present in P.
%
%   If supplied, this column is removed before fitting because glmnet
%   estimates its own intercept.
%
% useParallel
%   Fit neurons using parfor. Default: false.
%
% storeCVFits
%   Store all cvglmnet output objects. Default: false.
%
% maximumLinearPredictor
%   Upper numerical limit applied before exponentiation to prevent
%   floating-point overflow. Default: 700.
%
% OUTPUT
% ------
% model.K
%   E-by-N kernel matrix in log-rate space.
%
% model.intercept
%   1-by-N fitted log-rate intercepts.
%
% model.linearPredictor
%   T-by-N fitted log expected responses.
%
% model.Fhat
%   T-by-N fitted positive mean responses.
%
% model.deviance
%   Poisson deviance for each neuron.
%
% model.nullDeviance
%   Poisson deviance for an intercept-only model.
%
% model.devianceExplained
%   1 - deviance/nullDeviance for each neuron.
%
% model.rSquared
%   Response-scale squared-error R-squared. Retained as a descriptive
%   metric; deviance explained is usually preferable for Poisson models.
%
% model.lambda
%   Selected lambda for each neuron.
%
% model.cvFit
%   Cross-validation objects if storeCVFits is true.
%
% Notes
% -----
% Kernel coefficients may be negative. A negative coefficient represents
% a reduction in the log expected response. The final prediction remains
% positive because of the exponential link.
%
% Deconvolved activity may be continuous rather than integer-valued. In
% that case, the Poisson model is interpreted as a working likelihood or
% quasi-likelihood rather than a literal spike-count model.

arguments
    P (:,:) double
    F (:,:) double

    options.alpha (1,1) double = 0.01
    options.lambdaChoice (1,1) string = "lambda_min"
    options.nFolds (1,1) double {mustBeInteger, mustBePositive} = 5
    options.foldId (:,1) double = []

    options.standardizePredictors (1,1) logical = true
    options.includeGLMNetIntercept (1,1) logical = true
    options.interceptColumn (:,1) double = []

    options.useParallel (1,1) logical = false
    options.storeCVFits (1,1) logical = false

    options.maximumLinearPredictor (1,1) double = 700
end

%% Validate dimensions and options

[T, E] = size(P);
[responseT, N] = size(F);

if responseT ~= T
    error("fitDirectGLMNet:TimeMismatch", ...
        "P and F must have the same number of rows.");
end

if options.alpha < 0 || options.alpha > 1
    error("fitDirectGLMNet:InvalidAlpha", ...
        "alpha must lie between zero and one.");
end

if ~ismember(options.lambdaChoice, ["lambda_min", "lambda_1se"])
    error("fitDirectGLMNet:InvalidLambdaChoice", ...
        "lambdaChoice must be lambda_min or lambda_1se.");
end

if options.maximumLinearPredictor <= 0
    error("fitDirectGLMNet:InvalidMaximumLinearPredictor", ...
        "maximumLinearPredictor must be positive.");
end

if numel(options.interceptColumn) > 1
    error("fitDirectGLMNet:MultipleIntercepts", ...
        "interceptColumn must be empty or a single column index.");
end

if ~isempty(options.interceptColumn)
    if options.interceptColumn < 1 || ...
            options.interceptColumn > E || ...
            options.interceptColumn ~= round(options.interceptColumn)

        error("fitDirectGLMNet:InvalidInterceptColumn", ...
            "interceptColumn must index one column of P.");
    end
end

%% Validate Poisson responses

finiteResponses = F(isfinite(F));

if any(finiteResponses < 0)
    minimumResponse = min(finiteResponses);

    error("fitDirectGLMNet:NegativeResponse", ...
        ["Poisson regression requires nonnegative responses. " ...
         "The minimum observed response is %.6g."], ...
        minimumResponse);
end

%% Remove rows containing invalid values

validRows = all(isfinite(P), 2) & all(isfinite(F), 2);

if ~all(validRows)
    warning("fitDirectGLMNet:InvalidRows", ...
        "Removing %d rows containing NaN or Inf.", ...
        sum(~validRows));
end

Pfit = P(validRows,:);
Ffit = F(validRows,:);

Tfit = size(Pfit,1);

if Tfit < 3
    error("fitDirectGLMNet:TooFewObservations", ...
        "Too few valid observations remain.");
end

%% Remove an explicit intercept column

includedPredictorColumns = true(1,E);

if ~isempty(options.interceptColumn)
    includedPredictorColumns(options.interceptColumn) = false;
end

X = Pfit(:,includedPredictorColumns);
nPredictorsFit = size(X,2);

if isempty(X)
    error("fitDirectGLMNet:NoPredictors", ...
        "No predictors remain after removing the intercept column.");
end

%% Construct contiguous temporal folds

if isempty(options.foldId)
    foldId = makeContiguousFolds(Tfit, options.nFolds);
else
    suppliedFoldId = options.foldId(:);

    if numel(suppliedFoldId) == T
        foldId = suppliedFoldId(validRows);

    elseif numel(suppliedFoldId) == Tfit
        foldId = suppliedFoldId;

    else
        error("fitDirectGLMNet:FoldSizeMismatch", ...
            ["foldId must contain one value per original row or one " ...
             "value per retained valid row."]);
    end
end

if any(~isfinite(foldId)) || any(foldId < 1)
    error("fitDirectGLMNet:InvalidFoldId", ...
        "foldId must contain finite positive fold labels.");
end

% Convert arbitrary fold labels to consecutive integer labels.
[~, ~, foldId] = unique(foldId, "stable");

nFolds = numel(unique(foldId));

if nFolds < 3
    error("fitDirectGLMNet:TooFewFolds", ...
        "At least three cross-validation folds are required.");
end

%% Configure glmnet

glmOptions = glmnetSet;

glmOptions.alpha = options.alpha;
glmOptions.standardize = options.standardizePredictors;
glmOptions.intr = options.includeGLMNetIntercept;

%% Fit every neuron

Kfit = zeros(nPredictorsFit,N);
intercept = zeros(1,N);
selectedLambda = nan(1,N);

allZeroResponse = all(Ffit == 0,1);

if options.storeCVFits
    cvFit = cell(1,N);
else
    cvFit = {};
end

if options.useParallel
    if options.storeCVFits
        parfor neuronIndex = 1:N
            [Kfit(:,neuronIndex), ...
                intercept(neuronIndex), ...
                selectedLambda(neuronIndex), ...
                cvFit{neuronIndex}] = fitSingleNeuronDirectPoisson( ...
                    X, ...
                    Ffit(:,neuronIndex), ...
                    glmOptions, ...
                    nFolds, ...
                    foldId, ...
                    options.lambdaChoice);
        end
    else
        parfor neuronIndex = 1:N
            [Kfit(:,neuronIndex), ...
                intercept(neuronIndex), ...
                selectedLambda(neuronIndex)] = ...
                fitSingleNeuronDirectPoisson( ...
                    X, ...
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
            [Kfit(:,neuronIndex), ...
                intercept(neuronIndex), ...
                selectedLambda(neuronIndex), ...
                cvFit{neuronIndex}] = ...
                fitSingleNeuronDirectPoisson( ...
                    X, ...
                    Ffit(:,neuronIndex), ...
                    glmOptions, ...
                    nFolds, ...
                    foldId, ...
                    options.lambdaChoice);
        else
            [Kfit(:,neuronIndex), ...
                intercept(neuronIndex), ...
                selectedLambda(neuronIndex)] = ...
                fitSingleNeuronDirectPoisson( ...
                    X, ...
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

%% Restore original E-row kernel organization

K = zeros(E,N);
K(includedPredictorColumns,:) = Kfit;

% The removed intercept column has a kernel coefficient of zero because
% glmnet's separate intercept is used instead.
if ~isempty(options.interceptColumn)
    K(options.interceptColumn,:) = 0;
end

%% Generate log-link predictions

% This is eta = log(E[F | X]).
linearPredictorFit = X * Kfit + intercept;

% Guard against overflow in exp().
linearPredictorForExp = min( ...
    linearPredictorFit, ...
    options.maximumLinearPredictor);

FhatFit = exp(linearPredictorForExp);

% Represent all-zero neurons with a numerically tiny positive prediction.
FhatFit(:,allZeroResponse) = realmin;

%% Return predictions aligned to original rows

linearPredictor = nan(T,N);
linearPredictor(validRows,:) = linearPredictorFit;

Fhat = nan(T,N);
Fhat(validRows,:) = FhatFit;

%% Calculate Poisson deviance

deviance = poissonDeviance(Ffit, FhatFit);

% Intercept-only null model.
nullMean = mean(Ffit,1);
nullPrediction = repmat(nullMean, Tfit, 1);

% An all-zero neuron has a null mean of zero. Use realmin so that the
% numerical prediction remains strictly positive.
nullPrediction(:,nullMean == 0) = realmin;

nullDeviance = poissonDeviance(Ffit, nullPrediction);

devianceExplained = 1 - deviance ./ nullDeviance;
devianceExplained(nullDeviance == 0) = NaN;

populationDeviance = sum(deviance, "omitnan");
populationNullDeviance = sum(nullDeviance, "omitnan");

if populationNullDeviance > 0
    populationDevianceExplained = ...
        1 - populationDeviance / populationNullDeviance;
else
    populationDevianceExplained = NaN;
end

%% Calculate response-scale R-squared

residual = Ffit - FhatFit;

sumSquaredError = sum(residual.^2,1);

centeredResponse = Ffit - mean(Ffit,1);
totalSumSquares = sum(centeredResponse.^2,1);

rSquared = 1 - sumSquaredError ./ totalSumSquares;
rSquared(totalSumSquares == 0) = NaN;

populationDenominator = sum(centeredResponse(:).^2);

if populationDenominator > 0
    populationR2 = 1 - ...
        sum(residual(:).^2) / populationDenominator;
else
    populationR2 = NaN;
end

%% Calculate Poisson log likelihood without constant log(y!)

poissonLogLikelihood = sum( ...
    Ffit .* log(max(FhatFit,realmin)) - FhatFit, ...
    1);

%% Package outputs

model = struct;

model.family = "poisson";
model.link = "log";

model.K = K;
model.Kfit = Kfit;

model.intercept = intercept;

model.linearPredictor = linearPredictor;
model.Fhat = Fhat;

model.deviance = deviance;
model.nullDeviance = nullDeviance;
model.devianceExplained = devianceExplained;
model.populationDevianceExplained = ...
    populationDevianceExplained;

model.rSquared = rSquared;
model.populationR2 = populationR2;

model.poissonLogLikelihood = poissonLogLikelihood;

model.lambda = selectedLambda;
model.lambdaChoice = options.lambdaChoice;
model.alpha = options.alpha;

model.foldId = foldId;
model.validRows = validRows;

model.includedPredictorColumns = find(includedPredictorColumns);
model.removedInterceptColumn = options.interceptColumn;
model.allZeroResponse = allZeroResponse;

model.cvFit = cvFit;

end


function [weights, intercept, lambda, cvFit] = ...
    fitSingleNeuronDirectPoisson( ...
        X, response, glmOptions, nFolds, foldId, lambdaChoice)
%FITSINGLENEURONDIRECTPOISSON Fit one Poisson elastic-net model.
%
% The fitted model is:
%
%     log(E[y | X]) = intercept + X*weights

if any(~isfinite(response)) || any(response < 0)
    error("fitSingleNeuronDirectPoisson:InvalidResponse", ...
        "Poisson responses must be finite and nonnegative.");
end

% An identically zero neuron has no finite maximum-likelihood intercept.
% Return a tiny constant expected response instead.
if all(response == 0)
    weights = zeros(size(X,2),1);
    intercept = log(realmin);
    lambda = NaN;
    cvFit = [];
    return
end

cvFit = cvglmnet( ...
    X, ...
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

% glmnet returns the intercept in the first row.
intercept = full(coefficients(1));
weights = full(coefficients(2:end));

switch lambdaChoice
    case "lambda_min"
        lambda = cvFit.lambda_min;

    case "lambda_1se"
        lambda = cvFit.lambda_1se;
end

end


function deviance = poissonDeviance(observed, predicted)
%POISSONDEVIANCE Calculate column-wise Poisson deviance.
%
% For observed values greater than zero:
%
%   D = 2 * sum(y*log(y/mu) - y + mu)
%
% When y = 0, the y*log(y/mu) term is defined as zero.

if any(~isfinite(predicted(:))) || any(predicted(:) <= 0)
    error("poissonDeviance:InvalidPrediction", ...
        "Predictions must be finite and strictly positive.");
end

term = predicted - observed;

positiveObserved = observed > 0;

term(positiveObserved) = ...
    observed(positiveObserved) .* ...
    log(observed(positiveObserved) ./ ...
        predicted(positiveObserved)) ...
    - observed(positiveObserved) ...
    + predicted(positiveObserved);

deviance = 2 * sum(term,1);

end


function foldId = makeContiguousFolds(T, nFolds)
%MAKECONTIGUOUSFOLDS Split time points into contiguous temporal blocks.

if nFolds < 3
    error("makeContiguousFolds:TooFewFolds", ...
        "At least three folds are required.");
end

if nFolds > T
    error("makeContiguousFolds:TooManyFolds", ...
        "nFolds cannot exceed the number of observations.");
end

foldEdges = round(linspace(0,T,nFolds + 1));
foldId = zeros(T,1);

for foldIndex = 1:nFolds
    firstIndex = foldEdges(foldIndex) + 1;
    lastIndex = foldEdges(foldIndex + 1);

    foldId(firstIndex:lastIndex) = foldIndex;
end

end