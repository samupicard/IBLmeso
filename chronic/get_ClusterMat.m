function [outMatrix, outPs, outVs, outPos, clusterLookup, sessionLookup] = get_ClusterMat(dataStruct)

%from a ROI structure that contains the fields 'subject', 'date', 
%'session', 'clusterUID' and 'peth', produce a matrix of [nROIs,
%nConds, nSamples, nSessions] containing PETHs for each unique cluster
%and session.

% Filter out entries with empty clusterUIDs
isValid = arrayfun(@(d) ~isempty(d.clusterUID), dataStruct);
dataStruct = dataStruct(isValid);
K = numel(dataStruct);

% Extract identifiers
clusterUIDs = [dataStruct.clusterUID];
sessionKeys = arrayfun(@(d) sprintf('%s_%s_%s', d.subject, d.date, d.session), dataStruct, 'UniformOutput', false);

% Get unique clusters and sessions
[clusterLookup, ~, clusterIndices] = unique(clusterUIDs);        % [m x 1]
[sessionLookup, ~, sessionIndices] = unique(sessionKeys);        % [q x 1]
m = numel(clusterLookup);
n = numel(sessionLookup);

% Get size of PETH matrix
sz = size(dataStruct(1).peth);  % assuming all are same size and the last dimension is time

% Initialize output matrices
if length(sz)==3
    outMatrix = NaN(m, n, sz(3), sz(1), sz(2));
elseif length(sz)==2
    outMatrix = NaN(m, n, sz(2), sz(1));
end
outPs = NaN(m,n);

% Fill the matrix
for k = 1:K
    i = clusterIndices(k);  % index for cluster
    j = sessionIndices(k);  % index for session

    if ndims(outMatrix)==4
        outMatrix(i, j, :, :) = permute(dataStruct(k).peth, [2,1]);
    elseif ndims(outMatrix)==5
        outMatrix(i, j, :, :, :) = permute(dataStruct(k).peth, [3,1,2]);
    end

    % get corresponding p-values
    outPs(i,j) = dataStruct(k).p;
    
    % get corresponding value of the test-statistic
    outVs(i,j) = dataStruct(k).tstat_empirical;
    
    % get position for each cluster
    outPos(i,j,:) = dataStruct(k).pos;
    
end


end
