function [rShuffle, nearestIdx, nearestDist] = nearestNeighborRepeatReliability(tr1, tr2, fov, pos)
% For each ROI, correlate repeat 1 of that ROI with repeat 2 of nearest
% neighboring ROI in the same FOV.

nROIs = size(tr1,2);

rShuffle = nan(1,nROIs);
nearestIdx = nan(1,nROIs);
nearestDist = nan(1,nROIs);

fov = fov(:);
pos = pos(:,1:2);   % use XY only

for f = unique(fov(:))'

    idx = find(fov == f);

    if numel(idx) < 2
        continue
    end

    D = pdist2(pos(idx,:), pos(idx,:));
    D(1:numel(idx)+1:end) = inf;  % exclude self

    [nearestDistThis, iNearestLocal] = min(D, [], 2);

    nearestIdxThis = idx(iNearestLocal);

    nearestIdx(idx) = nearestIdxThis;
    nearestDist(idx) = nearestDistThis;

    rShuffle(idx) = columnwisePearson( ...
        tr1(:,idx), ...
        tr2(:,nearestIdxThis));
end

end