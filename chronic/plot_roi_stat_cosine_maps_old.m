function plot_roi_stat_cosine_maps_old(out, clim)
if nargin < 2 || isempty(clim)
    clim = [-1 1];
end

nSub = 5; %numel(out.subjects);
nLag = 16;%numel(out.lags);
nCol = nLag + 1;   % add within-day column

bas = aratopdown.atlas.build_topdown;


%figure('Color','w', 'Units','normalized', 'Position',[0.05 0.05 0.9 0.85]);
figure('Color','w','Position',[3000 50 1800 750]);
subCnt = 0;

for isub = 5+(1:nSub)

    subCnt = subCnt+1;

    % Precompute per-subject axis limits
    mask = false(size(out.meanCos{isub,1}));
    for ik = 1:nLag
        M = out.meanCos{isub,ik};
        mask = mask | ~isnan(M);
    end
    [rows, cols] = find(mask);
    if isempty(rows)
        xlims = [0,5200];
        ylims = [-4000,2000];
    else
        xlims = out.xEdges([min(cols),max(cols)]) + [-100,+100];
        ylims = out.yEdges([min(rows),max(rows)]) + [-100,+100];
    end

    for icol = 1:nCol

        if icol == 1
            M = out.withinDayMeanCos{isub};
        else
            ik = icol - 1;
            M = out.meanCos{isub,ik};
        end

        if any(~isnan(M(:))) && ~isempty(M)

            subplot(nSub, nLag, (subCnt-1)*nLag + ik);

            hImg = imagesc(out.xEdges(1:end-1), out.yEdges(1:end-1), M);
            axis image tight
            set(gca, 'YDir', 'normal');

            % NaNs → transparent
            set(hImg, 'AlphaData', ~isnan(M));

            % Background color for NaNs
            ax = gca;
            ax.Color = [0.7 0.7 0.7];

            caxis(clim);
            colormap(brewermap([],'*RdBu'))

            hold on
            cellfun(@(x) cellfun(@(x) plot(1000*x(:,2),1000*x(:,1),'color',[0 0 0]),x,'uni',false), ...
                {bas.dorsal_brain_areas(1:end-11).boundaries_stereotax},'uni', false);
            %xlim([0,5200]); ylim([-5100,3600]);
            %xlim([0,5200]); ylim([-4000,2000]);
            xlim(xlims); ylim(ylims);
            set(ax,'xtick',[],'ytick',[]);
            if ik > 1
                ax.YAxis.Visible = 'off';
            end
            ax.XAxis.Visible = 'off';
            axis square
            daspect([1 1 1])

            box off

        end

        if ik == 1
            ylabel(char(out.subjects(isub)));
        end
    end
end

cb = colorbar('Position',[0.92 0.11 0.015 0.77]);
cb.Label.String = 'Cosine similarity';
end