function cropped_ims = crop_cluster_ims(ims)
    % Crops the images to the smallest rectangle containing all non-zero pixels.
    %
    % Args:
    %     ims (3D array): Images to crop (size: [n, H, W])
    %
    % Returns:
    %     cropped_ims (3D array): Cropped images (size: [n, H', W'])

    % Compute the max projection over all images (collapse first dimension)
    ims_max = squeeze(max(ims, [], 1));

    % Find nonzero pixels
    [rows, cols] = find(ims_max > 0);

    % Determine cropping boundaries
    z_top = max(rows);
    z_bottom = min(rows);
    z_left = min(cols);
    z_right = max(cols);

    % Ensure boundaries stay within image limits
    z_bottom = max(z_bottom - 1, 1);
    z_top = min(z_top + 1, size(ims, 2));
    z_left = max(z_left - 1, 1);
    z_right = min(z_right + 1, size(ims, 3));

    % Crop the images while preserving the first dimension (n)
    cropped_ims = ims(:, z_bottom:z_top, z_left:z_right);

    % Set the border pixels to 1
    %cropped_ims(:, [1, end], :) = 1;
    %cropped_ims(:, :, [1, end]) = 1;
end
