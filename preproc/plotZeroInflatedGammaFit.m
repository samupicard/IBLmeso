function plotZeroInflatedGammaFit(x, nBins, col, varargin)

p = inputParser;
p.addParameter('UpperQuantile', 0.95, @isscalar);
p.addParameter('Epsilon', 1e-4, @isscalar);
p.parse(varargin{:});

x = x(:);
x = x(isfinite(x));
xPos = x(x > 0);

if isempty(xPos)
    return
end

zig = fitZeroInflatedGamma(x, 'Epsilon', p.Results.Epsilon);

xMaxPlot = quantile(xPos, p.Results.UpperQuantile);
edges = linspace(0, xMaxPlot, nBins + 1);

histogram(xPos, edges, ...
    'Normalization','probability', ...
    'EdgeColor','none', ...
    'FaceColor',col, ...
    'FaceAlpha',0.5);

hold on;

if isfinite(zig.k) && isfinite(zig.a)
    % Expected conditional probability per bin under Gamma
    binProb = gamcdf(edges(2:end) + zig.epsilon, zig.k, zig.a) - ...
              gamcdf(edges(1:end-1) + zig.epsilon, zig.k, zig.a);

    binCtr = edges(1:end-1) + diff(edges)/2;

    plot(binCtr, binProb, 'k', 'LineWidth', 2);
end

txt = sprintf('q = %.2f\nloc = %.3g\nk = %.2f\na = %.2g', ...
    zig.q, zig.loc, zig.k, zig.a);

text(0.98, 0.98, txt, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top');

xlabel('value');
ylabel('probability');
box off;

end