function [theta, R2] = computePairMetrics(x, y)
% compute theta + R^2 for one population-vector pair

theta = NaN;
R2 = NaN;

if numel(x) < 2 || numel(y) < 2
    return
end

% theta
denom = norm(x) * norm(y);
if denom > 0
    c = dot(x, y) / denom;
    c = max(min(c, 1), -1); % numerical safety
    theta = acosd(c);
end

% R^2 from linear regression y ~ x
pfit = polyfit(x, y, 1);
yhat = polyval(pfit, x);
ssRes = sum((y - yhat).^2);
ssTot = sum((y - mean(y)).^2);
if ssTot > 0
    R2 = 1 - ssRes/ssTot;
end
end


