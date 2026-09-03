function plot_roi_stat_fracSig_maps(out, clim, bas)
% Plot per-subject fraction-significant maps from roi_stat_fracSig_maps.
%
% Usage:
%   plot_roi_stat_fracSig_maps(out)
%   plot_roi_stat_fracSig_maps(out, [0 1])
%   plot_roi_stat_fracSig_maps(out, [0 0.5], bas)

if nargin < 2 || isempty(clim)
    clim = [0 1];
end

if nargin < 3
    bas = aratopdown.atlas.build_topdown;
end

nSub = 4;%numel(out.subjects);

% Choose a roughly square layout
nCol = ceil(sqrt(nSub));
nRow = ceil(nSub / nCol);

nCol = 1;
nRow = nSub;

figure('Color','w','Position',[3000 50 300 750]);
subCnt = 0;

for isub = [3,4,9,10]%1:nSub

    subCnt = subCnt + 1;

    % Per-subject axis limits
    mask = ~isnan(out.meanFracSig{isub});
    [rows, cols] = find(mask);

    if isempty(rows)
        xlims = [0, 5200];
        ylims = [-4000, 2000];
    else
        xlims = out.xEdges([min(cols), max(cols)]) + [-100, +100];
        ylims = out.yEdges([min(rows), max(rows)]) + [-100, +100];
    end
    xlims = [400, 4800];
    ylims = [-4100, 1400];

    M = out.meanFracSig{isub};

    subplot(nRow, nCol, subCnt);

    if any(~isnan(M(:))) && ~isempty(M)

        hImg = imagesc(out.xEdges(1:end-1), out.yEdges(1:end-1), M);
        axis image tight
        set(gca, 'YDir', 'normal');

        % NaNs -> transparent
        set(hImg, 'AlphaData', ~isnan(M));

        % Background color for NaNs
        ax = gca;
        ax.Color = [0.7 0.7 0.7];

        caxis(clim);
        %colormap(parula);
        colormap(brewermap([],'YlOrRd'));

        hold on
        if ~isempty(bas) && isfield(bas, 'dorsal_brain_areas')
            cellfun(@(x) cellfun(@(y) plot(1000*y(:,2),1000*y(:,1), ...
                'color', [0 0 0]), x, 'uni', false), ...
                {bas.dorsal_brain_areas(1:end-11).boundaries_stereotax}, ...
                'uni', false);
        end

        xlim(xlims);
        ylim(ylims);
        set(ax, 'xtick', [], 'ytick', []);
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';

        axis square
        daspect([1 1 1])
        box off

    else
        ax = gca;
        ax.Color = [0.7 0.7 0.7];
        xlim(xlims);
        ylim(ylims);
        set(ax, 'xtick', [], 'ytick', []);
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';
        axis square
        daspect([1 1 1])
        box off
    end

    title(char(out.subjects(isub)), 'Interpreter', 'none');
end

cb = colorbar('Position',[0.75 0.11 0.09 0.77]);
cb.Label.String = 'Fraction significant';