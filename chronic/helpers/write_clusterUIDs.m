function cUIDs_all = write_clusterUIDs(subject, dates)

%write_chronicUUIDs takes a ROICaT output file,
%assigns a UUID to each FOV / cluster combination, and then writes
%corresponding arrays of UUIDs to each of the alf collections
%
%written by Samuel Picard

root = 'Y:\Subjects\';
datpath = fullfile(root,subject,'Chronic');

save_csv = true; %set this to false if you just want to test

% find the ROICaT data
fns = dir(fullfile(datpath,'*ROICaT.tracking.results.mat')); %needs to contain mpci data

%TODO prepare an output named cUIDs_all, and fill this with the cUIDs of
%each cFOV.

for iFOV = 1:length(fns)
    fprintf('cFOV %d/%d : ', iFOV, length(fns))

    %load the ROICaT output data
    load(fullfile(fns(iFOV).folder,fns(iFOV).name));

    % Check whether sample silhouette metrics exist
    has_QMs = false;
    if exist('quality_metrics', 'var')
        if isfield(quality_metrics, 'sample_silhouette') && isfield(quality_metrics, 'cluster_silhouette')
            sample_silhouette_all = [quality_metrics.sample_silhouette{:}];
            cluster_silhouettes = [quality_metrics.cluster_silhouette{:}];
            has_QMs = true;
        end
    end

    %TODO make a selection of sessions based on dates
    % (make a cell array of date strings from input_data.paths_stat, and
    % compare this to the cell array of dates)
    selected_sessions = true(1,length(input_data.paths_stat));
    nsessions = sum(selected_sessions);

    selected_sessions_idxs = find(selected_sessions);

    %generate UUIDs
    nclusters = size(clusters.labels_bool,2)-1;
    cUIDs = strings(nclusters,1);
    for i = 1:nclusters
        cUIDs(i) = char(java.util.UUID.randomUUID);
    end

    %save these in the Chronic folder
    if save_csv
        %save(fullfile(datpath,'cluster_uids_all.mat'), 'uuids', '-v7.3');  % v7.3 supports large arrays
        writematrix(cUIDs,fullfile(datpath,[fns(iFOV).name(1:end-28),'.clusterUIDs_all.csv']));
    end

    %now assign cluster ID to each ROI in each session
    fprintf('Session ');
    nCharsPrinted = 0;
    sess_cntr = 0;
    roi_cntr = 0;

    for iSession = 1:nsessions

        fprintf(repmat('\b', 1, nCharsPrinted))
        nCharsPrinted = fprintf('%d/%d ..', iSession, nsessions);

        sessionIdx = selected_sessions_idxs(iSession);

        cidxs = clusters.labels_bySession{sessionIdx};
        nROIs = length(cidxs);

        if has_QMs
            sample_silhouettes = sample_silhouette_all(roi_cntr + (1:nROIs));
        end

        % make arrays by ROI
        CUIDs = strings(nROIs,1);
        if has_QMs
            ss_byROI = nan(nROIs,1);
            cs_byROI = nan(nROIs, 1);
        end

        for iROI = 1:nROIs
            if cidxs(iROI) > -1
                CUIDs(iROI) = cUIDs(cidxs(iROI)+1);
                if has_QMs
                    ss_byROI(iROI) = sample_silhouettes(iROI);
                    cs_byROI(iROI) = cluster_silhouettes(cidxs(iROI)+1);
                end
            end
        end

        roi_cntr = roi_cntr + nROIs;

        % save output to alf folder
        pth = fullfile(input_data.paths_stat{sessionIdx},'..','..');
        pth_new = replace(pth,'128.40.198.163','128.40.198.150');
        path_alf = dir(pth_new);
        ff = path_alf(3).folder;

        % assert that the arrays match the mpciROIs length
        ctx = sprintf('Subject=%s, FOV=%s, sessionIdx=%d', subject, fns(iFOV).name, sessionIdx);
        assert_cUIDs_matches_mpciROIs(ff, nROIs, ctx);

        if save_csv
            writematrix(CUIDs, fullfile(ff,'mpciROIs.clusterUIDs.csv'));
            if has_QMs
                writeNPY(ss_byROI, fullfile(ff,'mpciROIs.sampleSilhouette.npy'));
                writeNPY(cs_byROI, fullfile(ff,'mpciROIs.clusterSilhouette.npy'));
            end
        end
    end

    fprintf('. Done!\n');

end

end

function assert_cUIDs_matches_mpciROIs(alfFolder, expectedN, contextStr)
% Assert mpciROIs.mpciROITypes.npy exists in alfFolder and matches expectedN.
%
% expectedN: typically nROIs from clusters.labels_bySession{iSession}
% contextStr: string for nicer error messages (optional)

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
    'ROICaT:CUIDLengthMismatch: %s\nExpected nROIs=%d but mpciROIs.mpciROITypes has %d entries.\nFolder: %s', ...
    contextStr, expectedN, nCC, alfFolder);
end

