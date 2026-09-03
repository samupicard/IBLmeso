function [is_tracked, iROIs] = get_trackedROIs(clusterUIDs,sessionPaths)

%get_trackedROIs  returns indices of ROIs tracked across multiple sessions  
%
% get_trackedROIs(clusterUIDs,sessionPaths) takes an array of clusterUIDs 
% (strings), and returns a boolean array of the same size, which is true 
% for those UIDs which are present in all the sessions indicated in the 
% cell array sessionPaths. 
% 
% [is_tracked, iROIs] = get_trackedROIs(clusterUIDs,sessionPaths) also 
% returns an array of indices (sorted in alphabetical order of clusterUIDs).
%
% written by Samuel Picard

%load the UIDs for all other sessions
cluster_ids_all = {};
for iSess = 1:length(sessionPaths)
    
    datpath = sessionPaths{iSess};
    cluster_ids_all{iSess} = get_clusterUIDs(datpath);
        
end

%get non-empty UIDs
shared_uids = clusterUIDs(clusterUIDs ~= "");  % exclude empties

%find the intersection across all the other UIDs
for j = 1:numel(cluster_ids_all)
    shared_uids = intersect(shared_uids, cluster_ids_all{j});
end

% Return boolean mask
is_tracked = ismember(clusterUIDs, shared_uids);

% Return a set of indices, ordered alphabetically by clusterUIDs
[~, iROIs] = ismember(shared_uids, clusterUIDs);