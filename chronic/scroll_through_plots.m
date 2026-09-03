function scroll_through_plots()
    % Create sample data (100 timepoints x 50 plots)
    nPlots = 50;
    data = randn(100, nPlots); 

    % Set up the UIFigure
    fig = uifigure('Name', 'Plot Navigator', 'Position', [100 100 700 500]);

    % Set up axes
    ax = uiaxes(fig, ...
        'Position', [50 100 600 350], ...
        'XGrid', 'on', 'YGrid', 'on');
    xlabel(ax, 'Time');
    ylabel(ax, 'Value');

    % Set up Next and Previous buttons
    btnPrev = uibutton(fig, ...
        'Text', 'Previous', ...
        'Position', [200 30 100 30], ...
        'ButtonPushedFcn', @(btn,event) navigatePlot(-1));

    btnNext = uibutton(fig, ...
        'Text', 'Next', ...
        'Position', [400 30 100 30], ...
        'ButtonPushedFcn', @(btn,event) navigatePlot(1));

    % Initialize current index
    currentIndex = 1;
    plotCurrent();

    % Callback to update the plot
    function navigatePlot(direction)
        % Update index within bounds
        currentIndex = currentIndex + direction;
        currentIndex = max(1, min(nPlots, currentIndex));
        plotCurrent();
    end

    % Plot the current data column
    function plotCurrent()
        cla(ax);
        plot(ax, data(:, currentIndex), '-o');
        title(ax, sprintf('Plot %d of %d', currentIndex, nPlots));
    end
end
