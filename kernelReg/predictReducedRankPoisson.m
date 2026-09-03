function Fhat = predictReducedRankPoisson(model, Pnew)

    % Apply the training-set centering
    Pc = Pnew - model.predictorMean;

    % Poisson log-rate
    eta = Pc * model.K + model.intercept;

    % Prevent numerical overflow in exp()
    eta = min(eta, 700);

    % Predicted mean response
    Fhat = exp(eta);

end