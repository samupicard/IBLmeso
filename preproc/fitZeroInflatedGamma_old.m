function zig = fitZeroInflatedGamma_old(x)
%FITZEROINFLATEDGAMMA Fit zero-inflated shifted gamma model 
%
% Model:
%   s ~ (1-q)*delta(0) + q*Gamma(k,a,loc)
%
% Here:
%   q   = P(s > 0)
%   loc = minimum positive value
%   k   = gamma shape
%   a   = gamma scale
%
% written by Samuel Picard based on Wei,... & Paninski (2020)

x = x(:);
x = x(isfinite(x));

xPos = x(x > 0);

zig.nTotal = numel(x);
zig.nPos = numel(xPos);
zig.q = zig.nPos / zig.nTotal;

zig.loc = NaN;
zig.k = NaN;
zig.a = NaN;

if numel(xPos) < 10
    return
end

zig.loc = min(xPos);

y = xPos - zig.loc;
y = y(y > 0);

if numel(y) < 10
    return
end

mu = mean(y);
v  = var(y);

if v <= 0 || mu <= 0
    return
end

% Gamma method-of-moments estimates
zig.k = mu^2 / v;   % shape
zig.a = v / mu;     % scale

end