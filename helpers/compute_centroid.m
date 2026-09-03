function centers = compute_centroid(sparse_data)
    % COMPUTE_CENTROID Computes the center of mass for N 2D objects stored as sparse rows.
    %
    % Args:
    %     sparse_data (sparse double, N x 262144): Flattened sparse images.
    %
    % Returns:
    %     centers (N x 2 double): Center of mass coordinates [row, col] for each object.
    
    % Grid size
    grid_size = 512;
    N = size(sparse_data, 1); % Number of objects
    
    % Precompute grid indices (row, col) for all 512x512 locations
    [X, Y] = meshgrid(1:grid_size, 1:grid_size);
    
    % Flatten indices to match sparse format (1 x 262144)
    X = X(:)';
    Y = Y(:)';
    
    % Initialize output
    centers = zeros(N, 2);
    
    % Compute center of mass for each object
    for i = 1:N
        % Extract nonzero values and indices
        [cols, values] = find(sparse_data(i, :));
        
        % Compute weighted sum of coordinates
        total_mass = sum(values);
        if total_mass > 0
            centers(i, 1) = sum(values .* Y(cols)) / total_mass; % Row (y) coordinate
            centers(i, 2) = sum(values .* X(cols)) / total_mass; % Column (x) coordinate
        else
            centers(i, :) = NaN; % Handle empty objects
        end
    end
end
