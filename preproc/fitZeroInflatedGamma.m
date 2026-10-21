function zig = fitZeroInflatedGamma(x, varargin)
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

p = inputParser;
p.addParameter('Epsilon', 1e-4, @isscalar);
p.addParameter('MinNPositive', 10, @isscalar);
p.parse(varargin{:});

eps0 = p.Results.Epsilon;

x = x(:);
x = x(isfinite(x));

xPos = x(x > 0);

zig.nTotal = numel(x);
zig.nPos   = numel(xPos);
zig.q      = zig.nPos / zig.nTotal;

zig.loc = NaN;
zig.k   = NaN;
zig.a   = NaN;
zig.epsilon = eps0;

if numel(xPos) < p.Results.MinNPositive
    return
end

% Descriptive minimum positive value only
zig.loc = min(xPos);

% Gamma fit to positive values, nudged away from zero
y = xPos + eps0;

warnState = warning('off','MATLAB:nearlySingularMatrix');
pd = fitdist(y(:),'Gamma');
warning(warnState);

zig.k = pd.a;
zig.a = pd.b;


end