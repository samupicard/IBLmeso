function nearest_neighbours = find_nearest_neighbours(sparse_data)
    % FIND_NEAREST_NEIGHBOURS Finds the nearest object for each object based on centroid distance.
    %
    % Args:
    %     sparse_data (sparse double, N x [H*W]): Flattened sparse images.
    %
    % Returns:
    %     nearest_neighbors (N x 1 double): Index of the nearest neighbor for each object.

    % Compute centroids
    centers = compute_centroids(sparse_data); % [N x 2] matrix of [row, col] coordinates

    % Get number of objects
    N = size(centers, 1);

    % Initialize output
    nearest_neighbours = NaN(N, 1);

    % Compute pairwise distances and find nearest neighbors
    for i = 1:N
        % Compute Euclidean distance to all other centroids
        distances = sqrt(sum((centers - centers(i, :)).^2, 2));

        % Set self-distance to infinity so it's not counted
        distances(i) = inf;

        % Find the index of the minimum distance (nearest neighbor)
        [~, nearest_neighbours(i)] = min(distances);
    end
end
