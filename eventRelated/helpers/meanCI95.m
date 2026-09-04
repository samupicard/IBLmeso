function [mu, lo, hi] = meanCI95(vals)
% mean and 95% CI using normal approximation

vals = vals(isfinite(vals));

if isempty(vals)
    mu = NaN; lo = NaN; hi = NaN;
    return
end

mu = mean(vals);

if numel(vals) == 1
    lo = mu;
    hi = mu;
    return
end

sem = std(vals, 0, 1) / sqrt(numel(vals));
delta = 1.96 * sem;

lo = mu - delta;
hi = mu + delta;
end