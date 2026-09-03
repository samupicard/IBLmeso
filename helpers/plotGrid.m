function plotGrid(x, y, v)
% x, y: coordinates (N x 1)
% v: associated values (N x 1), between 0 and 1

    % Grid settings
    gridSize = 100; % spacing between grid lines
    min_x = floor(min(x)/gridSize)*gridSize;
    max_x = ceil(max(x)/gridSize)*gridSize;
    min_y = floor(min(y)/gridSize)*gridSize;
    max_y = ceil(max(y)/gridSize)*gridSize;
    edges_x = min_x:gridSize:max_x; 
    edges_y = min_y:gridSize:max_y;
    
    % Bin each point into x and y bins
    [~, xBin] = histc(x, edges_x);
    [~, yBin] = histc(y, edges_y);

    % Grid size (number of cells)
    numBins_x = length(edges_x) - 1;
    numBins_y = length(edges_y) - 1;

    % Initialize accumulators
    sumGrid = zeros(numBins_y, numBins_x);
    countGrid = zeros(numBins_y, numBins_x);

    % Accumulate p values in each bin
    for i = 1:length(v)
        xi = xBin(i);
        yi = yBin(i);
        if xi >= 1 && xi <= numBins_x && yi >= 1 && yi <= numBins_y
            sumGrid(yi, xi) = sumGrid(yi, xi) + v(i);     % row = y, col = x
            countGrid(yi, xi) = countGrid(yi, xi) + 1;
        end
    end

    % Compute average p in each bin
    avgGrid = sumGrid ./ countGrid;
    avgGrid(isnan(avgGrid)) = NaN;  % Handle empty bins

    % Plot the grid
    figure;
    imagesc(edges_x(1:end-1)+gridSize/2, edges_y(1:end-1)+gridSize/2, avgGrid);
    axis xy;
    %colorbar;
    %xlabel('x');
    %ylabel('y');
    %title('Average p value in each 5x5 grid cell');

end
