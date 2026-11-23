rootPath = 'Y:\Subjects';
jumpThreshold = 1e5;

% Recursively find all wheel position files
%files = dir(fullfile(rootPath, '**', '_ibl_wheel.position.npy'));

% Check each file for large jumps
pathsWithJumps = {};
progressMsg = '';

for i = 1:numel(files)
    fprintf(repmat('\b', 1, numel(progressMsg)));
    progressMsg = sprintf('Checking file %d/%d (%.1f%%)', ...
        i, numel(files), 100*i/numel(files));
    fprintf('%s', progressMsg);

    filePath = fullfile(files(i).folder, files(i).name);

    try
        wheelPos = readNPY(filePath);

        if any(abs(diff(double(wheelPos(:)))) >= jumpThreshold)
            pathsWithJumps{end+1,1} = filePath;
        end

    catch ME
        warning('Could not read %s: %s', filePath, ME.message);
    end
end

fprintf('Found %d files with jumps >= %g (out of %d files).\n', ...
    numel(pathsWithJumps), jumpThreshold, numel(files));

sessionPaths = regexprep(pathsWithJumps, ...
    '^.*?Subjects\\([^\\]+)\\([^\\]+)\\([^\\]+)\\.*$', '$1/$2/$3');

disp(sessionPaths);