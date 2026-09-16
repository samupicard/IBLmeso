
H = 512;
W = 512;
nFramesToLoad = 25000;
rootPath = 'Y:\Subjects\SP081';
fn_bin = 'imaging.frames_motionRegistered.bin';

%% find registered binary paths

% Find date folders
dateDirs = dir(rootPath);
dateDirs = dateDirs([dateDirs.isdir]);
%

% Keep yyyy-mm-dd folders only
isDate = ~cellfun('isempty', ...
    regexp({dateDirs.name}, '^\d{4}-\d{2}-\d{2}$', 'once'));

dateDirs = dateDirs(isDate);

dateDirs = dateDirs(end,:); %last day

files = [];

for iDate = 1:numel(dateDirs)

    datePath = fullfile(dateDirs(iDate).folder, dateDirs(iDate).name);

    % Find session folders, e.g. 001, 002
    sessDirs = dir(datePath);
    sessDirs = sessDirs([sessDirs.isdir]);

    isSession = ~cellfun('isempty', ...
        regexp({sessDirs.name}, '^\d{3}$', 'once'));

    sessDirs = sessDirs(isSession);

    for iSess = 1:numel(sessDirs)

        suite2pPath = fullfile( ...
            sessDirs(iSess).folder, ...
            sessDirs(iSess).name, ...
            'suite2p');

        if ~isfolder(suite2pPath)
            continue
        end

        % Find plane folders
        planeDirs = dir(fullfile(suite2pPath, 'plane*'));
        planeDirs = planeDirs([planeDirs.isdir]);

        for iPlane = 1:numel(planeDirs)

            f = dir(fullfile( ...
                planeDirs(iPlane).folder, ...
                planeDirs(iPlane).name, ...
                fn_bin));

            if ~isempty(f)
                files = [files; f]; %#ok<AGROW>
            end
        end
    end
end

fprintf('Found %d files.\n', numel(files));

%% Preallocate results

nFiles = numel(files);

filePath           = strings(nFiles, 1);
saturatedRatio     = nan(nFiles, 1);
globalMax          = nan(nFiles, 1);
maxTimeSaturated   = nan(nFiles, 1);

%% Run QC

for i = 1:nFiles

    filePath(i) = fullfile(files(i).folder, files(i).name);

    fprintf('[%d/%d] %s\n', i, nFiles, filePath(i));

    try
        [saturatedRatio(i), maxImg, timeSaturatedImg] = ...
            getSaturatedPixelRatio( ...
                filePath(i), H, W, nFramesToLoad);

        globalMax(i) = max(maxImg, [], 'all');
        maxTimeSaturated(i) = max(timeSaturatedImg, [], 'all');

    catch ME
        warning('Failed to process %s:\n%s', ...
            filePath(i), ME.message);
    end

end

%% Collect results

results = table( ...
    filePath, ...
    saturatedRatio, ...
    globalMax, ...
    maxTimeSaturated);

%% Display

disp(results);