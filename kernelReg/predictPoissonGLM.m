function Fpred = predictPoissonGLM(model, P)
%PREDICTPOISSONGLM Compute predicted Poisson mean activity.
%
% Output:
%   Fpred is time × neurons.

if isfield(model, "linearPredictor") && ...
        size(model.linearPredictor, 1) == size(P, 1)

    % Safest option when predicting the original fitted dataset.
    eta = model.linearPredictor;

elseif isfield(model, "rawPredictorIntercept")
    % Uncentered representation.
    eta = P * model.K + model.rawPredictorIntercept;

elseif isfield(model, "predictorMean")
    % Centered representation.
    Pc = P - model.predictorMean;
    eta = Pc * model.K + model.intercept;

else
    % Assumes P and model.K are already aligned.
    eta = P * model.K + model.intercept;
end

% Avoid numerical overflow.
eta = min(eta, 700);

Fpred = exp(eta);

end