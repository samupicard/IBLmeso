function roi_stat_cosine_gui(out, bas, mapField, clim, subjectSel)
% Interactive explorer for ROI cosine maps.
%
% Click any map bin to open/update a lag-timecourse plot for that
% subject and spatial bin.
%
% Inputs:
%   out      : output from roi_stat_cosine_maps
%   bas         : atlas/boundary struct (optional)
%   mapField    : "meanCos", "binShuffleP", "dayShuffleP"
%   clim        : color limits for displayed maps
%   subjectSel  : optional subject selection:
%                   - [] or omitted: all subjects
%                   - numeric indices, e.g. [2 4 5]
%                   - cellstr/string array, e.g. {'SP058','SP060'}
%
% Example:
%   roi_stat_cosine_gui(out_st, bas, "meanCos", [0 1], {'SP058','SP060'})

if nargin < 2
    bas = [];
end
if nargin < 3 || isempty(mapField)
    mapField = "meanCos";
end
if nargin < 4 || isempty(clim)
    clim = [0 1];
end
if nargin < 5
    subjectSel = [];
end

mapField = string(mapField);

allSubjects = string(out.subjects(:));

% Resolve subject selection
if isempty(subjectSel)
    subjIdx = 1:numel(allSubjects);
elseif isnumeric(subjectSel)
    subjIdx = subjectSel(:)';
elseif iscell(subjectSel) || isstring(subjectSel)
    wanted = string(subjectSel(:));
    subjIdx = find(ismember(allSubjects, wanted));
else
    error('subjectSel must be empty, numeric indices, or subject names.');
end

if isempty(subjIdx)
    error('No matching subjects found.');
end

showSubjects = allSubjects(subjIdx);
nSub = numel(subjIdx);
nLag = numel(out.lags);
lags = out.lags;

figMaps = figure( ...
    'Color','w', ...
    'Position',[3000 50 1800 150], ...
    'Name','ROI cosine maps', ...
    'NumberTitle','off');

axList = gobjects(nSub*nLag,1);
axIdx = 0;

for iSubPlot = 1:nSub
    isub = subjIdx(iSubPlot);

    % Subject-specific bounds from empirical maps
    mask = false(size(out.meanCos{isub,1}));
    for ik = 1:nLag
        Mtmp = out.meanCos{isub,ik};
        mask = mask | ~isnan(Mtmp);
    end
    [rows, cols] = find(mask);
    if isempty(rows)
        xlims = [400, 4800];
        ylims = [-4100, 1400];
    else
        xlims = out.xEdges([min(cols),max(cols)]) + [-300,+300];
        ylims = out.yEdges([min(rows),max(rows)]) + [-200,+800];
    end

    for icol = 1:nLag
        ik = icol;

        switch mapField
            case "meanCos"
                M = out.meanCos{isub,ik};
            case "binShuffleP"
                M = out.binShuffleP{isub,ik};
            case "dayShuffleP"
                M = out.dayShuffleP{isub,ik};
            otherwise
                error('Unsupported mapField: %s', mapField);
        end

        ax = subplot(nSub, nLag, (iSubPlot-1)*nLag + icol, 'Parent', figMaps);
        axIdx = axIdx + 1;
        axList(axIdx) = ax;

        if ~isempty(M) && any(~isnan(M(:)))
            [X,Y] = meshgrid(out.xEdges, out.yEdges);
            hImg = pcolor(X, Y, [M nan(size(M,1),1); nan(1,size(M,2)+1)]);
            shading flat
            set(ax, 'YDir', 'normal');
            %set(hImg, 'AlphaData', ~isnan(M));
            set(hImg, 'PickableParts', 'all');
            set(hImg, 'ButtonDownFcn', @(src,evt)clickMap(src, evt, ax, isub, out, figMaps));

            ax.Color = [0.7 0.7 0.7];
            caxis(ax, clim);
            colormap(ax, brewermap([],'*RdBu'));

            hold(ax, 'on');
            if ~isempty(bas) && isfield(bas, 'dorsal_brain_areas')
                cellfun(@(x) cellfun(@(y) plot(ax, 1000*y(:,2), 1000*y(:,1), ...
                    'Color',[0 0 0]), ywrap(x), 'uni', false), ...
                    {bas.dorsal_brain_areas(1:end-11).boundaries_stereotax}, ...
                    'uni', false);
            end

            xlim(ax, xlims);
            ylim(ax, ylims);
            set(ax,'XTick',[],'YTick',[]);
            if icol > 1
                ax.YAxis.Visible = 'off';
            end
            ax.XAxis.Visible = 'off';
            axis(ax, 'square');
            daspect(ax, [1 1 1]);
            box(ax, 'off');
        else
            ax.Color = [0.7 0.7 0.7];
            xlim(ax, xlims);
            ylim(ax, ylims);
            set(ax,'XTick',[],'YTick',[]);
            if icol > 1
                ax.YAxis.Visible = 'off';
            end
            ax.XAxis.Visible = 'off';
            axis(ax, 'square');
            daspect(ax, [1 1 1]);
            box(ax, 'off');
        end

        if iSubPlot == 1
            title(ax, sprintf('k = %d', lags(ik)));
        end
        if icol == 1
            ylabel(ax, char(out.subjects(isub)));
        end
    end
end

axList = axList(1:axIdx);
for i = 1:numel(axList)
    caxis(axList(i), clim);
end

cb = colorbar(axList(1), 'Position',[0.92 0.11 0.015 0.77]);
cb.Limits = clim;
cb.Ticks = linspace(clim(1), clim(2), 5);

switch mapField
    case "meanCos"
        cb.Label.String = 'Cosine similarity';
    case "binShuffleP"
        cb.Label.String = 'Bin-shuffle p-value';
    case "dayShuffleP"
        cb.Label.String = 'Day-shuffle p-value';
end

setappdata(figMaps, 'out_st', out);
setappdata(figMaps, 'selectedRect', gobjects(0));
setappdata(figMaps, 'selectedAxes', gobjects(0));
end

function clickMap(~, ~, ax, isub, out_st, figMaps)
cp = ax.CurrentPoint;
xClick = cp(1,1);
yClick = cp(1,2);

col = discretize(xClick, out_st.xEdges);
row = discretize(yClick, out_st.yEdges);

if isnan(row) || isnan(col)
    return
end

% Highlight selected bin on the clicked axes
oldRect = getappdata(figMaps, 'selectedRect');
oldAx   = getappdata(figMaps, 'selectedAxes');

if ~isempty(oldRect) && isgraphics(oldRect)
    delete(oldRect);
end
if ~isempty(oldAx) && isgraphics(oldAx)
    hold(oldAx, 'on');
end

x0 = out_st.xEdges(col);
x1 = out_st.xEdges(col+1);
y0 = out_st.yEdges(row);
y1 = out_st.yEdges(row+1);

hold(ax, 'on');
hRect = rectangle(ax, ...
    'Position', [x0, y0, x1-x0, y1-y0], ...
    'EdgeColor', [1 1 0], ...
    'LineWidth', 1);

setappdata(figMaps, 'selectedRect', hRect);
setappdata(figMaps, 'selectedAxes', ax);

nLag = numel(out_st.lags);
lags = out_st.lags;

vals = nan(1, nLag);
vals_dayShuff = nan(out_st.nShuff, nLag);
vals_binShuff = nan(out_st.nShuff, nLag);

for ik = 1:nLag
    M = out_st.meanCos{isub,ik};
    if ~isempty(M) && row <= size(M,1) && col <= size(M,2)
        vals(ik) = M(row,col);
    end

    D = out_st.dayShuffleCos{isub,ik};
    if ~isempty(D) && row <= size(D,1) && col <= size(D,2)
        vals_dayShuff(:,ik) = squeeze(D(row,col,:));
    end

    B = out_st.binShuffleCos{isub,ik};
    if ~isempty(B) && row <= size(B,1) && col <= size(B,2)
        vals_binShuff(:,ik) = squeeze(B(row,col,:));
    end
end

qs_dayShuff = quantile(vals_dayShuff,[0.025,0.5,0.975],1);
qs_binShuff = quantile(vals_binShuff,[0.025,0.5,0.975],1);

figTC = findobj('Type','figure', 'Name','ROI lag timecourse');
if isempty(figTC) || ~isvalid(figTC)
    figTC = figure('Color','w', 'Position',[3200 650 350 250], ...
        'Name','ROI lag timecourse', 'NumberTitle','off');
else
    figure(figTC);
    clf(figTC);
end

hold on
patch([lags, fliplr(lags)], ...
      [qs_dayShuff(1,:), fliplr(qs_dayShuff(end,:))], ...
      [.5 .5 .5], 'EdgeColor','none', 'FaceAlpha',0.25);
%plot(lags, qs_dayShuff(2,:), '--', 'Color',[.5 .5 .5], 'LineWidth', 2);
yline(0,'k--');
xline(0,'Color',[.5 .5 .5]);
% patch([lags, fliplr(lags)], ...
%       [qs_binShuff(1,:), fliplr(qs_binShuff(end,:))], ...
%       'r', 'EdgeColor','none', 'FaceAlpha',0.25);

plot(lags, vals, 'k-o', 'LineWidth', 2, 'MarkerSize', 4);

xlabel('k');
ylabel('Cosine similarity');
title(sprintf('%s | ML=[%.0f, %.0f] AP=[%.0f, %.0f]', ...
    char(out_st.subjects(isub)), ...
    out_st.xEdges(col), out_st.xEdges(col+1), ...
    out_st.yEdges(row), out_st.yEdges(row+1)));

%legend({'day shuffle 95% CI', 'bin shuffle 95% CI', 'empirical'}, ...
%    'Location', 'best');
grid off
box off

yAll = [vals(:); qs_dayShuff(:); qs_binShuff(:)];
yAll = yAll(isfinite(yAll));
if ~isempty(yAll)
    yPad = 0.05 * max(1e-6, range(yAll));
    ylim([min(yAll)-yPad, max(yAll)+yPad]);
end
ylim([-0.4,0.85]);
xlim([lags(1),lags(end)]);
end

function out = ywrap(x)
if iscell(x)
    out = x;
else
    out = {x};
end
end