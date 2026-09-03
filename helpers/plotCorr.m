function h = plotCorr(xx,yy,ii)

if nargin<3
    ii=false(1,length(xx));
end

[R,P] = corrcoef(xx,yy);

cols = colororder;

hold on;
h(1) = scatter(xx,yy,'MarkerEdgeColor',cols(1,:));
h(2) = scatter(xx(ii),yy(ii),'MarkerEdgeColor',cols(1,:),'MarkerFaceColor',cols(1,:));
ylims = get(gca,'Ylim');
xlims = get(gca,'Xlim');

h(3) = text(xlims(1)+0.9*diff(xlims),ylims(1)+0.8*diff(ylims),sprintf('R%c=%.2f\np=%.3f',178,R(2)^2,P(2)),...
    'FontSize',10,'HorizontalAlignment','right','Color','w');