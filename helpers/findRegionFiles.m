function [regions, paths] = findRegionFiles(folder, wildcardPattern)
%findRegionFiles finds files matching wildcardPattern*.mat and extracts REGION from _REGION.mat
%
% Example wildcardPattern: 'GoodSessions_x_y_z_*'
% This function will search for: fullfile(folder, [wildcardPattern '*.mat'])

    files = dir(fullfile(folder, [wildcardPattern '*.mat']));

    regions = {};
    paths   = {};

    for k = 1:numel(files)
        name = files(k).name;

        % REGION is everything after the last underscore, before .mat
        tok = regexp(name, '_([^_]+)\.mat$', 'tokens');
        if isempty(tok)
            % doesn't follow the expected suffix naming rule, ignore
            continue
        end

        regions{end+1} = tok{1}{1}; %#ok<AGROW>
        paths{end+1}   = fullfile(folder, name); %#ok<AGROW>
    end

    % Guard against duplicates (same region appears multiple times in same folder/pattern)
    if numel(unique(regions)) ~= numel(regions)
        [u, ~, ic] = unique(regions);
        counts = accumarray(ic, 1);
        dup = u(counts > 1);
        error('Duplicate region(s) found for pattern "%s": %s', wildcardPattern, strjoin(dup, ', '));
    end
end