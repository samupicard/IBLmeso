function slopes = get_R2_decay_slopes(R2s_byType, varargin)
% get_R2_decay_slopes
%
% Fit a linear slope of R^2 vs lag (k) for each pethType.
%
% Input:
%   R2s_byType : cell array, one entry per pethType
%       each entry must contain:
%           .mean   (1 x K) vector of mean R^2 vs lag k=1..K
%
% Optional:
%   'minPoints' : minimum number of valid points required (default = 2)
%
% Output:
%   slopes : [nTypes x 1] slope per pethType

p = inputParser;
p.addParameter('minPoints', 2, @(x)isnumeric(x) && isscalar(x) && x>=2);
p.parse(varargin{:});
opt = p.Results;

nTypes = numel(R2s_byType);
slopes = nan(nTypes,1);

for s = 1:nTypes
    if isempty(R2s_byType{s}) || ~isfield(R2s_byType{s}, 'mean')
        continue
    end

    y = R2s_byType{s}.mean(:);   % column
    K = numel(y);
    x = (1:K)';                  % lag values

    keep = isfinite(x) & isfinite(y);

    if sum(keep) < opt.minPoints
        continue
    end
    
    %take out last point from fit
    x_nolast = x(keep); x_nolast = x_nolast(1:end-1);
    y_nolast = y(keep); y_nolast = y_nolast(1:end-1);

    pfit = polyfit(x_nolast, y_nolast, 1);
    slopes(s) = pfit(1);   % slope
end
end