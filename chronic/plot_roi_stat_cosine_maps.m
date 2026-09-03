function plot_roi_stat_cosine_maps(out, plotWhat, clims, bas)
% Plot ROI cosine similarity maps from roi_stat_cosine_maps output.
%
% Supports:
%   "meanCos"
%   "binShuffleMean"
%   "binShuffleP"
%   "binShuffleZ"
%   "dayShuffleMean"
%   "dayShuffleP"
%   "dayShuffleZ"
%
% Usage:
%   plot_roi_stat_cosine_maps(out)
%   plot_roi_stat_cosine_maps(out, "meanCos", [-1 1], bas)
%   plot_roi_stat_cosine_maps(out, "binShuffleP", [0 0.2], bas)
%   plot_roi_stat_cosine_maps(out, "dayShuffleZ", [0 5], bas)

if nargin < 2 || isempty(plotWhat)
    plotWhat = "meanCos";
end
plotWhat = string(plotWhat);

if nargin < 3 || isempty(clims)
    switch plotWhat
        case {"meanCos","binShuffleMean","dayShuffleMean"}
            clims = [-1 1];
        case {"binShuffleP","dayShuffleP"}
            clims = [0 1];
        case {"binShuffleZ","dayShuffleZ"}
            clims = [0 5];
        otherwise
            error('Unknown plotWhat: %s', plotWhat);
    end
end

if nargin < 4
    bas = aratopdown.atlas.build_topdown;
end

nSub = 1;%numel(out.subjects);
nLag = numel(out.lags);
nCol = nLag;

figure('Color','w','Position',[3000 50 1800 150]);
subCnt = 0;
axList = gobjects(nSub*nCol,1);
axIdx = 0;

for isub = 4%[3,4,9,10]%1:nSub
    subCnt = subCnt + 1;

    % Per-subject bounds from empirical maps
    mask = false(size(out.meanCos{isub,1}));
    for ik = 1:nLag
        Mtmp = out.meanCos{isub,ik};
        mask = mask | ~isnan(Mtmp);
    end
    [rows, cols] = find(mask);
    if isempty(rows)
        xlims = [0,5200];
        ylims = [-4000,2000];
    else
        xlims = out.xEdges([min(cols),max(cols)]) + [-100,+100];
        ylims = out.yEdges([min(rows),max(rows)]) + [-100,+100];
    end

    colcnt = 0;
    for icol = 1:nCol
        ik = icol;
        colcnt = colcnt+1;

        switch plotWhat
            case "meanCos"
                M = out.meanCos{isub,ik};
            case "binShuffleMean"
                M = out.binShuffleMean{isub,ik};
            case "binShuffleP"
                M = out.binShuffleP{isub,ik};
            case "binShuffleZ"
                M = out.binShuffleZ{isub,ik};
            case "dayShuffleMean"
                M = out.dayShuffleMean{isub,ik};
            case "dayShuffleP"
                M = out.dayShuffleP{isub,ik};
            case "dayShuffleZ"
                M = out.dayShuffleZ{isub,ik};
            otherwise
                error('Unknown plotWhat: %s', plotWhat);
        end

        subplot(nSub, nCol, (subCnt-1)*nCol + icol);
        axIdx = axIdx + 1;
        axList(axIdx) = gca;

        if ~isempty(M) && any(~isnan(M(:)))
            hImg = imagesc(out.xEdges(1:end-1), out.yEdges(1:end-1), M);
            axis image tight
            set(gca, 'YDir', 'normal');

            % NaNs -> transparent
            set(hImg, 'AlphaData', ~isnan(M));

            % Background color for NaNs
            ax = gca;
            ax.Color = [0.7 0.7 0.7];

            clim(clims);
            switch plotWhat
                case {"meanCos","binShuffleMean","dayShuffleMean"}
                    colormap(brewermap([],'*RdYlBu'));
                case {"binShuffleP","dayShuffleP"}
                    colormap(brewermap([],'*Reds'));
                otherwise
                    colormap(brewermap([],'*YlOrRd'));
            end

            hold on
            if ~isempty(bas) && isfield(bas, 'dorsal_brain_areas')
                cellfun(@(x) cellfun(@(y) plot(1000*y(:,2),1000*y(:,1), ...
                    'color',[0 0 0]), ywrap(x), 'uni', false), ...
                    {bas.dorsal_brain_areas(1:end-11).boundaries_stereotax}, ...
                    'uni', false);
            end

            xlim(xlims);
            ylim(ylims);
            set(ax,'xtick',[],'ytick',[]);

            if icol > 1
                ax.YAxis.Visible = 'off';
            end
            ax.XAxis.Visible = 'off';

            axis square
            daspect([1 1 1])
            box off
        else
            ax = gca;
            ax.Color = [0.7 0.7 0.7];
            xlim(xlims);
            ylim(ylims);
            set(ax,'xtick',[],'ytick',[]);
            if icol > 1
                ax.YAxis.Visible = 'off';
            end
            ax.XAxis.Visible = 'off';
            axis square
            daspect([1 1 1])
            box off
        end

        if subCnt == 1
            title(sprintf('k = %d', out.lags(ik)));
        end

        if icol == 1
            ylabel(char(out.subjects(isub)));
        end
    end
end

for i = 1:axIdx
    clim(axList(i), clims);
end

cb = colorbar(axList(1),'Position',[0.92 0.11 0.015 0.77]);
switch plotWhat
    case "meanCos"
        cb.Label.String = 'Cosine similarity';
    case "binShuffleMean"
        cb.Label.String = 'Bin-shuffle cosine';
    case "binShuffleP"
        cb.Label.String = 'Bin-shuffle p-value';
    case "binShuffleZ"
        cb.Label.String = 'Bin-shuffle Z';
    case "dayShuffleMean"
        cb.Label.String = 'Day-shuffle cosine';
    case "dayShuffleP"
        cb.Label.String = 'Day-shuffle p-value';
    case "dayShuffleZ"
        cb.Label.String = 'Day-shuffle Z';
end
end

function out = ywrap(x)
if iscell(x)
    out = x;
else
    out = {x};
end
end
           
