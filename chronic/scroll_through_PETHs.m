function scroll_through_PETHs(data,Cs,Ss,T)

    % Parameters
    nClusters = size(data,1);       % Number of pages (each page has 6 plots)
    nConds = size(data,2);
    nSamples = size(data,3);  % Length of each individual plot
    if ndims(data)==4
        nSessions = size(data,4);
    else
        nSessions = 1;
    end
    
    % Create synthetic data (a 3D array: [time x plot x page])
    data = randn(nClusters, nSamples, nConds, nSessions); 
    
    

    % Set up UIFigure
    fig = uifigure('Name', 'PETH navigator', 'Position', [100 100 1200 400]);

    % Set up axes in a row
    axesArray = gobjects(1, nSessions);
    for i = 1:nSessions
        axesArray(i) = uiaxes(fig, ...
            'Position', [50 + (i-1)*180, 80, 160, 250], ...
            'XTickLabelMode', 'auto', ...
            'YTickLabelMode', 'auto', ...
            'XGrid', 'on', ...
            'YGrid', 'on');
    end

    % Buttons
    btnPrev = uibutton(fig, ...
        'Text', 'Previous', ...
        'Position', [400 20 100 30], ...
        'ButtonPushedFcn', @(btn,event) navigatePage(-1));

    btnNext = uibutton(fig, ...
        'Text', 'Next', ...
        'Position', [520 20 100 30], ...
        'ButtonPushedFcn', @(btn,event) navigatePage(1));

    % Initialize current page index
    currentROI = 1;
    plotCurrentPage();

    % Navigation callback
    function navigatePage(direction)
        currentROI = currentROI + direction;
        currentROI = max(1, min(nClusters, currentROI));
        plotCurrentPage();
    end

    % Plot current page
    function plotCurrentPage()
        for i = 1:nSessions
            ax = axesArray(i);
            cla(ax);
            plot(ax, data(currentROI, :, i));
            title(ax, sprintf('Plot %d-%d', currentROI, i));
        end
    end
end

