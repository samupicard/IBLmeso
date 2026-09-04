function ax = IBL_plotROISummaryMap(T, varargin)

p = inputParser;
addParameter(p, 'alpha2', 0.025, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'sizeScale', 0.02, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'minSize', 1, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'binSize', 100, @(x)isnumeric(x) && isscalar(x)); % microns
addParameter(p, 'ax', [], @(x) isempty(x) || isgraphics(x,'axes'));
addParameter(p, 'bas', aratopdown.atlas.build_topdown, @isstruct);
addParameter(p, 'plotNonsigAsDots', false, @islogical);
addParameter(p, 'xLimits', [0 5200], @(x)isnumeric(x) && numel(x)==2);
addParameter(p, 'yLimits', [-5100 3600], @(x)isnumeric(x) && numel(x)==2);
parse(p, varargin{:});
opts = p.Results;

if isempty(opts.ax)
    figure('Color','k');
    ax = axes('Color','k');
else
    ax = opts.ax;
end
hold(ax, 'on');

ML = T.ML;
AP = T.AP;
pval = T.p;
stat = T.stat;

valid = isfinite(ML) & isfinite(AP) & isfinite(pval) & isfinite(stat);
ML = ML(valid);
AP = AP(valid);
pval = pval(valid);
stat = stat(valid);

%define shades of red, blue and gray
cBlue = [0.35 0.60 0.95];
cRed  = [0.95 0.45 0.45];
cGray = [0.55 0.55 0.55];

% ----- Fixed 100 x 100 micron density bins -----
xLimits = opts.xLimits;
yLimits = opts.yLimits;

xEdges = xLimits(1):opts.binSize:xLimits(2);
yEdges = yLimits(1):opts.binSize:yLimits(2);

% Ensure edges cover the full range exactly, even if not divisible by binSize
if xEdges(end) < xLimits(2)
    xEdges(end+1) = xLimits(2);
end
if yEdges(end) < yLimits(2)
    yEdges(end+1) = yLimits(2);
end

% Compute dot size & alpha
dotSize = opts.minSize + opts.sizeScale .* abs(stat);
dotAlpha = 1 - (1-abs(1.8*(pval-0.5))).^0.4;
dotAlpha = 0.8 * ones(size(pval)); %force equal alpha for each dot

% non-significant ROIs
isGray = pval > opts.alpha2 & pval < 1 - opts.alpha2;

if sum(isGray)>0
    if opts.plotNonsigAsDots

        scatter(ax, ML(isGray), AP(isGray), dotSize(isGray), ...
            cGray, 'filled', ...
            'MarkerFaceAlpha', 'flat', ...
            'MarkerEdgeAlpha', 0.1, ...
            'AlphaData',dotAlpha(isGray));

    else

        N = histcounts2(ML(isGray), AP(isGray), xEdges, yEdges);

        % Log transform so dense regions do not dominate
        Nlog = log10(N + 1);

        % Optional robust clipping of extreme density values
        hi = prctile(Nlog(Nlog > 0), 99);
        if isempty(hi) || hi == 0
            hi = 1;
        end
        Nlog(Nlog > hi) = hi;

        % Plot density map
        xCenters = xEdges(1:end-1) + diff(xEdges)/2;
        yCenters = yEdges(1:end-1) + diff(yEdges)/2;

        imagesc(ax, xCenters, yCenters, Nlog');
        set(ax, 'YDir', 'normal');

        % Low-saturation black-to-gray colormap
        cmap = gray(256);
        cmap = cmap * 0.65;   % max is gray, not white
        colormap(ax, cmap);
        clim(ax, [0 hi]);

    end
end
% Plot atlas boundaries
cellfun(@(x) cellfun(@(y) ...
    plot(ax, 1000*y(:,2), 1000*y(:,1), 'Color', [.6 .6 .6]), ...
    x, 'uni', false), ...
    {opts.bas.dorsal_brain_areas(1:end-11).boundaries_stereotax}, ...
    'uni', false);

% ----- Significant ROIs -----

isBlue = pval <= opts.alpha2;
isRed  = pval >= 1 - opts.alpha2;

scatter(ax, ML(isBlue), AP(isBlue), dotSize(isBlue), ...
    cBlue, 'filled', ...
    'MarkerFaceAlpha', 'flat', ...
    'MarkerEdgeAlpha', 0.1, ...
    'AlphaData',dotAlpha(isBlue));


scatter(ax, ML(isRed), AP(isRed), dotSize(isRed), ...
    cRed, 'filled', ...
    'MarkerFaceAlpha', 'flat', ...
    'MarkerEdgeAlpha', 0.1, ...
    'AlphaData',dotAlpha(isRed));

xlim(ax, xLimits);
ylim(ax, yLimits);

axis(ax, 'square');
daspect(ax, [1 1 1]);

set(ax, ...
    'XTick', [], ...
    'YTick', [], ...
    'Color', 'k', ...
    'Box', 'on', ...
    'LineWidth', 0.5, ...
    'XColor', [1 1 1], ...
    'YColor', [1 1 1]);

ax.TickLength = [0 0];

end