function centroids = compute_binary_centroids(sparse_data)
    % COMPUTE_BINARY_CENTROIDS Computes centroids of binary objects on a 512x512 grid.
    % 
    % Args:
    %     sparse_data (sparse double, N x 262144): Flattened sparse images.
    %
    % Returns:
    %     centroids (N x 2 double): [row, col] centroids for each object.

    % Get number of objects
    N = size(sparse_data, 1);

    % Initialize output
    centroids = NaN(N, 2);

    % Loop through each object
    for i = 1:N
        % Extract nonzero indices (convert to binary mask)
        [~,cols] = find(sparse_data(i, :)); 
        
        if isempty(cols)
            continue;  % Skip empty objects
        end

        % Convert linear indices to (row, col) on a 512x512 grid
        Y = mod(cols - 1, 512) + 1;  % Row indices
        X = floor((cols - 1) / 512) + 1;  % Column indices

        % Compute centroid as the mean of binary pixel locations
        centroids(i, :) = [mean(Y), mean(X)];
    end
end
