function deviance = poissonDeviance(observed, predicted)
%POISSONDEVIANCE Calculate column-wise Poisson deviance.
%
% For y > 0:
%
%     D = 2 * sum(y*log(y/mu) - (y-mu))
%
% For y = 0, the y*log(y/mu) term is defined as zero.

if any(predicted(:) <= 0) || any(~isfinite(predicted(:)))
    error("poissonDeviance:InvalidPrediction", ...
        "Poisson predictions must be finite and strictly positive.");
end

term = predicted - observed;

positiveObserved = observed > 0;

term(positiveObserved) = ...
    observed(positiveObserved) .* ...
    log(observed(positiveObserved) ./ predicted(positiveObserved)) ...
    - observed(positiveObserved) ...
    + predicted(positiveObserved);

deviance = 2 * sum(term, 1);

end
