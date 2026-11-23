function h = plotPMfit(psychoM, varargin)
% plotPMfit Plot psychometric data and fitted curves by block
%
% h = plotPMfit(psychoM)
% h = plotPMfit(psychoM,'Name',Value,...)
%
% Required fields in psychoM:
%   blockProbLeft   : 1 x nBlocks vector
%   contrastDiffs   : 1 x nContrasts vector
%   choiceTypes     : 1 x 2 vector
%   numTrials       : nBlocks x nContrasts
%   mean            : nBlocks x nChoices x nContrasts
%   bias            : scalar fitted intercept
%   sensitivity     : scalar fitted contrast coefficient
%   exponent        : scalar exponent q
%   deltabias       : scalar fitted block coefficient (optional if 1 block)
%   stat.covb       : covariance of GLM coefficients [b, v, db] (or [b, v])
%
% Notes:
% - This plots the empirical fraction of "right" choices, assumed to be
%   choiceTypes(2), consistent with getPsycho_fit.
% - Confidence intervals are approximate 95% Wald intervals on the linear
%   predictor, transformed through the logistic function.
%
% Optional name-value pairs:
%   'ax'            : target axes handle (default: new axes)
%   'colors'        : nBlocks x 3 colormap (default: lines / custom 2-block)
%   'lineWidth'     : width of fitted lines (default: 2.5)
%   'alphaShade'    : CI shading alpha (default: 0.18)
%   'markerEdgeColor': edge color for points (default: 'none')
%   'markerFaceAlpha': point alpha (default: 0.95)
%   'minMarkerSize' : minimum marker area (default: 20)
%   'maxMarkerSize' : maximum marker area (default: 140)
%   'xGrid'         : x grid for smooth curves (default: linspace(min,max,300))
%   'showLegend'    : true/false (default: true)
%   'xlabel'        : x label (default: 'Contrast difference')
%   'ylabel'        : y label (default: 'P(choose right)')
%   'title'         : axes title (default: '')
%
% Samuel Picard

% ----------------------------
% Parse inputs
% ----------------------------
ip = inputParser;
ip.addParameter('ax', [], @(x) isempty(x) || isgraphics(x,'axes'));
ip.addParameter('colors', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==3));
ip.addParameter('lineWidth', 2.5, @isscalar);
ip.addParameter('alphaShade', 0.18, @isscalar);
ip.addParameter('markerEdgeColor', 'none');
ip.addParameter('markerFaceAlpha', 0.95, @isscalar);
ip.addParameter('minMarkerSize', 20, @isscalar);
ip.addParameter('maxMarkerSize', 140, @isscalar);
ip.addParameter('xGrid', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('showLegend', true, @islogical);
ip.addParameter('xlabel', 'Contrast difference', @ischar);
ip.addParameter('ylabel', 'P(choose right)', @ischar);
ip.addParameter('title', '', @ischar);
ip.parse(varargin{:});
opt = ip.Results;

if isempty(opt.ax)
    figure;
    ax = axes;
else
    ax = opt.ax;
end
hold(ax,'on');

% ----------------------------
% Extract data
% ----------------------------

% Exclude 0.5 blocks from plotting
keepBlocks = abs(psychoM.blockProbLeft(:) - 0.5) > 1e-6;

probas  = psychoM.blockProbLeft(keepBlocks)';
xData   = psychoM.contrastDiffs(:)';
nBlocks = numel(probas);
nX       = numel(xData);

assert(numel(psychoM.choiceTypes)==2, ...
    'This function expects exactly 2 choice types.');

% empirical P(choose right): second choice type
yData = reshape(psychoM.mean(keepBlocks,2,:),nBlocks,[]); % nBlocks x nContrasts
nData = psychoM.numTrials(keepBlocks,:); % nBlocks x nContrasts


if size(yData,1) ~= nBlocks || size(yData,2) ~= nX
    error('psychoM.mean / psychoM.numTrials dimensions do not match blockProbLeft and contrastDiffs.');
end

% ----------------------------
% Colors
% ----------------------------
if isempty(opt.colors)
    if nBlocks == 2
        % close to your example
        colors = [0.72 0.42 0.05;   % warm brown/orange
                  0.42 0.39 0.60];  % muted purple
    else
        colors = lines(nBlocks);
    end
else
    colors = opt.colors;
    if size(colors,1) < nBlocks
        error('Provided colors must have at least nBlocks rows.');
    end
end

% ----------------------------
% Marker sizes scaled by trial count
% ----------------------------
nFlat = nData(:);
if all(isnan(nFlat)) || all(nFlat==nFlat(1))
    ms = repmat(mean([opt.minMarkerSize opt.maxMarkerSize]), size(nData));
else
    nMin = min(nFlat);
    nMax = max(nFlat);
    ms = opt.minMarkerSize + (nData - nMin) ./ max(eps, (nMax - nMin)) ...
        * (opt.maxMarkerSize - opt.minMarkerSize);
end

% ----------------------------
% Smooth x-grid
% ----------------------------
if isempty(opt.xGrid)
    xGrid = linspace(min(xData), max(xData), 300);
else
    xGrid = opt.xGrid(:)';
end

% ----------------------------
% Set up fitted model quantities
% ----------------------------
hasBiasBlocks = isfield(psychoM,'deltabias') && ~isnan(psychoM.deltabias) && nBlocks > 1;

b  = psychoM.bias;
v  = psychoM.sensitivity;
q  = psychoM.exponent;

if hasBiasBlocks
    db = psychoM.deltabias;
else
    db = 0;
end

% Coefficient covariance for CI
covb = [];
if isfield(psychoM,'stat') && isstruct(psychoM.stat) && isfield(psychoM.stat,'covb')
    covb = psychoM.stat.covb;
end

% transformed contrast
fXgrid = sign(xGrid) .* abs(xGrid).^q;

% Decide mapping from plotted blocks to B values in the fitted model
% Convention used here:
%   low pLeft  -> B = +1
%   high pLeft -> B = -1
%   0.5 block  -> B = 0
Bplot = zeros(1,nBlocks);
if hasBiasBlocks
    if nBlocks == 2
        [~,ord] = sort(probas,'ascend');
        Bplot(ord(1)) =  1;  % lower pLeft => more rightward bias
        Bplot(ord(2)) = -1;  % higher pLeft => more leftward bias
    else
        % any 0.5 block gets 0, lower-than-0.5 gets +1, higher-than-0.5 gets -1
        tol = 1e-9;
        Bplot(probas < 0.5 - tol) =  1;
        Bplot(probas > 0.5 + tol) = -1;
        Bplot(abs(probas - 0.5) <= tol) = 0;
    end
end

% ----------------------------
% Plot each block
% ----------------------------
h = struct;
h.ax = ax;
h.scatter = gobjects(nBlocks,1);
h.line = gobjects(nBlocks,1);
h.patch = gobjects(nBlocks,1);

for iBlock = 1:nBlocks

    thisColor = colors(iBlock,:);
    thisB = Bplot(iBlock);

    % --- fitted mean prediction
    eta = b + v .* fXgrid + db .* thisB;
    pHat = 1 ./ (1 + exp(-eta));

    % --- confidence interval (if covariance available)
    pLo = nan(size(pHat));
    pHi = nan(size(pHat));

    if ~isempty(covb)
        if hasBiasBlocks
            betaRow = [ones(numel(xGrid),1), fXgrid(:), thisB*ones(numel(xGrid),1)];
        else
            betaRow = [ones(numel(xGrid),1), fXgrid(:)];
        end

        % keep only the matching block of covb in case glmfit returned extras
        covbUse = covb(1:size(betaRow,2), 1:size(betaRow,2));

        seEta = sqrt(sum((betaRow * covbUse) .* betaRow, 2));
        etaLo = eta(:) - 1.96 * seEta;
        etaHi = eta(:) + 1.96 * seEta;

        pLo = 1 ./ (1 + exp(-etaLo));
        pHi = 1 ./ (1 + exp(-etaHi));

        % shaded CI
        xx = [xGrid, fliplr(xGrid)];
        yy = [pLo(:)', fliplr(pHi(:)')];
        h.patch(iBlock) = fill(ax, xx, yy, thisColor, ...
            'FaceAlpha', opt.alphaShade, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end

    % fitted line
    h.line(iBlock) = plot(ax, xGrid, pHat, ...
        'Color', thisColor, ...
        'LineWidth', opt.lineWidth);

    % empirical points
    xThis = xData;
    yThis = yData(iBlock,:);
    sThis = ms(iBlock,:);

    valid = ~isnan(xThis) & ~isnan(yThis) & ~isnan(sThis);
    h.scatter(iBlock) = scatter(ax, xThis(valid), yThis(valid), sThis(valid), ...
        'MarkerFaceColor', thisColor, ...
        'MarkerEdgeColor', opt.markerEdgeColor, ...
        'MarkerFaceAlpha', opt.markerFaceAlpha);
end

% ----------------------------
% Cosmetics
% ----------------------------
ylim(ax, [0 1]);
xlim(ax, [min(xGrid), max(xGrid)]);
xlabel(ax, opt.xlabel);
ylabel(ax, opt.ylabel);
title(ax, opt.title);
box(ax,'off');
set(ax,'TickDir','out');
grid(ax,'on');

% Legend
if opt.showLegend
    legStr = strings(nBlocks,1);
    for iBlock = 1:nBlocks
        legStr(iBlock) = sprintf('pLeft = %.1f', probas(iBlock));
    end
    legend(ax, h.line, legStr, 'Location', 'best', 'Box', 'off');
end

end