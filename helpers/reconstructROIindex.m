function [fov_id, iROI] = reconstructROIindex(alfDir)
% reconstructROIindex
%
% Reconstructs, for each row in allNeurons.PETH.npy, the originating FOV
% and the ROI index within that FOV.
%
% Inputs
%   alfDir : path to top-level alf directory (e.g. './alf')
%
% Outputs
%   fov_id : [N x 1] FOV number (e.g. 0,1,2,...)
%   iROI   : [N x 1] ROI index within FOV (python 0-based)
%
% Requires:
%   readNPY.m
%   https://github.com/cortex-lab/npy_matlab

    arguments
        alfDir (1,:) char
    end

    % Find all FOV directories
    fovDirsAll = dir(fullfile(alfDir, 'FOV_*'));
    fovDirsAll = fovDirsAll([fovDirsAll.isdir]);

    if isempty(fovDirsAll)
        error('No FOV_* directories found in %s', alfDir);
    end

    % Keep only FOVs that actually contain ROI type file
    fovDirs = [];
    fovNums = [];

    for i = 1:numel(fovDirsAll)
        fovName = fovDirsAll(i).name;

        tok = regexp(fovName, '^FOV_(\d+)$', 'tokens', 'once');
        if isempty(tok)
            continue
        end
        fovNum = str2double(tok{1});

        roiTypeFile = fullfile(alfDir, fovName, 'mpciROIs.mpciROITypes.npy');
        if exist(roiTypeFile, 'file')
            fovDirs = [fovDirs; fovDirsAll(i)];
            fovNums = [fovNums; fovNum];
        end
    end

    if isempty(fovDirs)
        error('No FOV folders with mpciROIs.mpciROITypes.npy found in %s', alfDir);
    end

    % Sort valid FOVs by numeric ID
    [fovNums, order] = sort(fovNums);
    fovDirs = fovDirs(order);

    % Preallocate as empty (cheap compared to PETH size)
    fov_id = [];
    iROI   = [];

    % Loop only over valid FOVs
    for i = 1:numel(fovDirs)
        fovName = fovDirs(i).name;
        fovNum  = fovNums(i);

        roiTypeFile = fullfile(alfDir, fovName, 'mpciROIs.mpciROITypes.npy');
        roiTypes = readNPY(roiTypeFile);
        roiTypes = roiTypes(:);  % ensure column

        % ROIs classified as neurons
        neuronRoiIdx = find(roiTypes == 1);  % MATLAB 1-based

        % Append mapping in concatenation order
        fov_id = [fov_id; repmat(fovNum, numel(neuronRoiIdx), 1)];
        iROI   = [iROI;   neuronRoiIdx - 1]; % python 0-based
    end

end
