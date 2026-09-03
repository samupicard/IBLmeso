function [bestMatchIdx, bestMatchCorr, C] = findBestColumnMatches(X, Y, varargin)
% Find, for each column of Y, the best matching column of X.
%
% Optional constraints:
%   'XFov', 'YFov' : ROI FOV vectors
%   'XPos', 'YPos' : nROI x nDim centroid matrices
%   'MaxDist'      : max centroid distance, e.g. 10

p = inputParser;
p.addParameter('XFov', [], @isnumeric);
p.addParameter('YFov', [], @isnumeric);
p.addParameter('XPos', [], @isnumeric);
p.addParameter('YPos', [], @isnumeric);
p.addParameter('MaxDist', inf, @isscalar);
p.parse(varargin{:});
opts = p.Results;

if size(X,1) ~= size(Y,1)
    error('X and Y must have the same number of rows.');
end

nX = size(X,2);
nY = size(Y,2);

% Demean and normalize columns
X = X - mean(X,1,'omitnan');
Y = Y - mean(Y,1,'omitnan');

Xnorm = sqrt(sum(X.^2,1,'omitnan'));
Ynorm = sqrt(sum(Y.^2,1,'omitnan'));

X = X ./ Xnorm;
Y = Y ./ Ynorm;

% Full correlation matrix: nX x nY
C = X' * Y;

C(Xnorm == 0,:) = NaN;
C(:,Ynorm == 0) = NaN;

% Build validity mask
validMatch = true(nX,nY);

if ~isempty(opts.XFov) && ~isempty(opts.YFov)
    validMatch = validMatch & (opts.XFov(:) == opts.YFov(:)');
end

if ~isempty(opts.XPos) && ~isempty(opts.YPos) && isfinite(opts.MaxDist)
    validDist = false(nX,nY);

    fovs = unique(opts.YFov(:))';

    for f = fovs
        ix = find(opts.XFov(:) == f);
        iy = find(opts.YFov(:) == f);

        if isempty(ix) || isempty(iy)
            continue
        end

        D = pdist2(opts.XPos(ix,:), opts.YPos(iy,:));
        validDist(ix,iy) = D <= opts.MaxDist;
    end

    validMatch = validMatch & validDist;
end

% Exclude invalid pairs
C(~validMatch) = NaN;

[bestMatchCorr,bestMatchIdx] = max(C,[],1,'omitnan');

% Handle ROIs with no valid match
noMatch = all(isnan(C),1);
bestMatchIdx(noMatch) = NaN;
bestMatchCorr(noMatch) = NaN;
end