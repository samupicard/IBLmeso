function plotGammaFitNonzero(x, nBins, Color)

xpos = x(x > 0);
xpos = xpos(isfinite(xpos));

if numel(xpos) < 10
    return
end

% Fit gamma
pd = fitdist(xpos(:),'Gamma');

% Histogram
histogram(xpos, nBins, ...
    'Normalization','pdf', ...
    'EdgeColor','none','FaceColor',Color);

hold on;

% Fitted PDF
xx = linspace(min(xpos), max(xpos), 500);
plot(xx, pdf(pd,xx), ...
    'Color', Color, ...
    'LineWidth', 2);

% Annotate parameters
txt = sprintf('k=%.2f\ntheta=%.2f', ...
    pd.a, pd.b);

xl = xlim;
yl = ylim;

text(xl(2), yl(2), txt, ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top');

end