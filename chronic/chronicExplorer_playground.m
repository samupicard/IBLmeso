f = uifigure('Name', 'Plot Navigator');
ax = uiaxes('Parent', f, 'Position', [50 80 500 300]);
btnNext = uibutton(f, 'Text', 'Next', 'Position', [450 20 100 30]);
btnPrev = uibutton(f, 'Text', 'Previous', 'Position', [330 20 100 30]);

currentIdx = 1;
nPlots = 50;
data = randn(100, nPlots);

plotAt(ax,currentIdx);

% Button callbacks
btnNext.ButtonPushedFcn = @(btn, event) navigate(1);
btnPrev.ButtonPushedFcn = @(btn, event) navigate(-1);

% Update function
function plotAt(ax,idx)
    cla(ax);
    plot(ax, data(:, idx));
    title(ax, sprintf('Plot #%d', idx));
end

function navigate(step)
    currentIdx = max(1, min(nPlots, currentIdx + step));
    plotAt(ax,currentIdx);
end
