function pethTypes = discoverPethTypesAcrossSession(fovDirs)
%DISCOVERPETHTYPESACROSSSESSION Discover task and passive PETH types.
%
% Recognizes files such as:
%
%   mpciROIs.PETHavgNorm_stimOn_contrastDiff.npy
%   mpciROIs.PETHavgNorm_stimOn_contrastDiff_odd.npy
%   mpciROIs.PETHavgNorm_valveOn_passive.npy
%
% Returns type names such as:
%
%   stimOn_contrastDiff
%   valveOn_passive

pethTypes = {};

prefix = 'mpciROIs.PETHavgNorm_';
suffix = '.npy';

for iFov = 1:numel(fovDirs)

    fovFolder = fullfile( ...
        fovDirs(iFov).folder, ...
        fovDirs(iFov).name);

    files = dir(fullfile( ...
        fovFolder, ...
        'mpciROIs.PETHavgNorm_*.npy'));

    for iFile = 1:numel(files)

        fileName = files(iFile).name;

        % Ignore time-vector files, should any share this prefix.
        if contains(fileName, '.timeValues')
            continue
        end

        typeName = erase(fileName, prefix);
        typeName = erase(typeName, suffix);

        % Collapse odd/even versions onto the underlying PETH type.
        typeName = regexprep(typeName, '_(odd|even)$', '');

        if ~isempty(typeName)
            pethTypes{end+1} = typeName; %#ok<AGROW>
        end
    end
end

pethTypes = unique(pethTypes, 'stable');

end