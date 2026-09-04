function IBL_plotDensityPanel(D, frac, frac_clim)

make_hsv_map = @(frac, fraclim, counts, vlim) ...
    hsv2rgb(cat(3, ...
        0.70 * (1 - min(max((frac - fraclim(1)) ./ diff(fraclim), 0), 1)), ...
        ones(size(frac)), ...
        min(max((counts - vlim(1)) ./ diff(vlim), 0), 1) .^ D.v_gamma));

rgb = make_hsv_map(frac, frac_clim, D.n_sess_perBin, D.v_clim_sess);
rgb(repmat(isnan(frac),1,1,3)) = 0;

image(D.xedges, D.yedges, permute(rgb,[2 1 3]));
set(gca, 'YDir', 'normal');
axis image
hold on

cellfun(@(x) cellfun(@(y) ...
    plot(1000*y(:,2),1000*y(:,1),'color',[.6 .6 .6]), ...
    x,'uni',false), ...
    {D.bas.dorsal_brain_areas(1:end-11).boundaries_stereotax}, ...
    'uni', false);

%full right hemisphere
xlim([0,5200]);
ylim([-5100,3600]);

%manual
%xlim([500,4000]);
%ylim([-1000,2500]);

axis square
set(gca, 'XTick',[], 'YTick',[], 'XColor','none', 'YColor','none', 'Color','k', 'Box','off');
daspect([1 1 1])

end