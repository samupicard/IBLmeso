function interactiveClusterViewer(clusterPETHs, T, cUIDs)
    % Inputs:
    % - clusterPETHs: 5D matrix [clusters, days, time, conditions, contrasts]
    % - T: time vector
    % - example_cUIDs: cell array of cluster labels

    % Setup
    nClusters = size(clusterPETHs, 1);
    nDays = size(clusterPETHs, 2);
    nContrasts = size(clusterPETHs, 5);
    nConditions = size(clusterPETHs,4);
    
    interval_toplot = 1:15;
    xlims = [T(interval_toplot(1)), T(interval_toplot(end))];
    colors = brewermap(nConditions, 'PuOr') * 0.8;

    currentIndex = 1;

    % Create figure
    f = figure('Name', 'Interactive Cluster Viewer', 'Position', [2, 210, 1000, 601]);

    % UI: Previous Button
    uicontrol('Style', 'pushbutton', ...
        'String', '< Previous', ...
        'Units', 'normalized', ...
        'Position', [0.01 0.01 0.12 0.05], ...
        'FontSize', 12, ...
        'Callback', @previousCallback);

    % UI: Next Button
    uicontrol('Style', 'pushbutton', ...
        'String', 'Next >', ...
        'Units', 'normalized', ...
        'Position', [0.87 0.01 0.12 0.05], ...
        'FontSize', 12, ...
        'Callback', @nextCallback);

    % UI: Cluster label
    label = uicontrol('Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.15 0.01 0.7 0.05], ...
        'String', '', ...
        'FontSize', 12, ...
        'HorizontalAlignment', 'center');

    % Layout for plots
    axLayout = tiledlayout(f, nContrasts, nDays, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    % Initial plot
    updatePlot();

    % --- Callback Functions ---
    function previousCallback(~, ~)
        if currentIndex > 1
            currentIndex = currentIndex - 1;
            updatePlot();
        end
    end

    function nextCallback(~, ~)
        if currentIndex < nClusters
            currentIndex = currentIndex + 1;
            updatePlot();
        end
    end

    % --- Update Plot Function ---
    function updatePlot()
        delete(axLayout.Children);  % Clear old plots
        label.String = sprintf('Cluster %d: %s', currentIndex, cUIDs(currentIndex));

        maxval = max(clusterPETHs(currentIndex,:,interval_toplot,:,:), [], 'all');

        for iC = 1:nContrasts
            for iD = 1:nDays
                ax = nexttile(axLayout);
                hold on;
                for iCond = 1:nConditions
                    plot(T, squeeze(clusterPETHs(currentIndex,iD,:,iCond,iC)), 'LineWidth', 2, 'Color', colors(iCond,:));
                end
                xline(0);
                xlim(xlims);
                ylim([0, maxval * 1.05]);
                box off;

                if iC ~= nContrasts
                    ax.XTick = [];
                end
                if iD ~= 1
                    ax.YTick = [];
                end
            end
        end

        %title(axLayout, cUIDs(currentIndex));
        xlabel(axLayout, 'Session');
        ylabel(axLayout, 'Contrast diff');
    end
end
