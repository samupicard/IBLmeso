function patch_clusterSampleSilhouette(subject, dates)

root = 'Y:\Subjects\';
datpath = fullfile(root, subject, 'Chronic');

save_npy = true;

fns = dir(fullfile(datpath, '*ROICaT.tracking.results.mat'));

for iFOV = 1:length(fns)
    fprintf('cFOV %d/%d : ', iFOV, length(fns))

    load(fullfile(fns(iFOV).folder, fns(iFOV).name));

    % Skip old ROICaT outputs without quality metrics
    has_QMs = false;
    if exist('quality_metrics', 'var')
        if isfield(quality_metrics, 'sample_silhouette') && isfield(quality_metrics, 'cluster_silhouette')
            sample_silhouette_all = [quality_metrics.sample_silhouette{:}];
            cluster_silhouettes = [quality_metrics.cluster_silhouette{:}];
            has_QMs = true;
        end
    end

    if ~has_QMs
        fprintf('No cluster / sample silhouette metrics found. Skipping.\n');
        continue
    end

    selected_sessions = true(1, length(input_data.paths_stat));
    nsessions = sum(selected_sessions);
    selected_sessions_idxs = find(selected_sessions);

    fprintf('Session ');
    nCharsPrinted = 0;
    roi_cntr = 0;

    for iSession = 1:nsessions

        fprintf(repmat('\b', 1, nCharsPrinted))
        nCharsPrinted = fprintf('%d/%d ..', iSession, nsessions);

        sessionIdx = selected_sessions_idxs(iSession);

        cidxs = clusters.labels_bySession{sessionIdx};
        nROIs = length(cidxs);

        sample_silhouettes = sample_silhouette_all(roi_cntr + (1:nROIs));        

        ss_byROI = nan(nROIs, 1);
        cs_byROI = nan(nROIs, 1);
        for iROI = 1:nROIs
            if cidxs(iROI) > -1
                ss_byROI(iROI) = sample_silhouettes(iROI);
                cs_byROI(iROI) = cluster_silhouettes(cidxs(iROI)+1);
            end
        end

        roi_cntr = roi_cntr + nROIs;

        % Find matching ALF folder
        pth = fullfile(input_data.paths_stat{sessionIdx}, '..', '..');
        pth_new = replace(pth, '128.40.198.163', '128.40.198.150');
        path_alf = dir(pth_new);
        ff = path_alf(3).folder;

        ctx = sprintf('Subject=%s, FOV=%s, sessionIdx=%d', ...
            subject, fns(iFOV).name, sessionIdx);

        assert_cUIDs_matches_mpciROIs(ff, nROIs, ctx);

        if save_npy

            % remove old incorrectly named file if it exists
            old_fn = fullfile(ff, 'mpciROIs.clusterSampleSilhouette.npy');
            if isfile(old_fn)
                delete(old_fn);
            end

            % write corrected filenames
            writeNPY(ss_byROI, fullfile(ff, 'mpciROIs.sampleSilhouette.npy'));
            writeNPY(cs_byROI, fullfile(ff, 'mpciROIs.clusterSilhouette.npy'));
        end

    end

    fprintf('. Done!\n');
end

end

function assert_cUIDs_matches_mpciROIs(alfFolder, expectedN, contextStr)

if nargin < 3, contextStr = ''; end

ccPath = fullfile(alfFolder, 'mpciROIs.mpciROITypes.npy');
if ~isfile(ccPath)
    error('ROICaT:MissingMpciROITypes', ...
        'Missing mpciROIs.mpciROITypes.npy in %s. %s', alfFolder, contextStr);
end

cc = readNPY(ccPath);
nCC = size(cc, 1);
if isvector(cc)
    nCC = numel(cc);
end

assert(nCC == expectedN, ...
    ['ROICaT:CUIDLengthMismatch: %s\n' ...
     'Expected nROIs=%d but mpciROIs.mpciROITypes has %d entries.\n' ...
     'Folder: %s'], ...
    contextStr, expectedN, nCC, alfFolder);
end