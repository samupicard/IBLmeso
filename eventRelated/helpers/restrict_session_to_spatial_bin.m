function Sbin = restrict_session_to_spatial_bin(S, mlRange, apRange)
% restrict_session_to_spatial_bin
%
% Keep only neurons whose [ML, AP] coordinates fall inside the bin.
%
% mlRange = [ml0 ml1]
% apRange = [ap0 ap1]

if ~isfield(S, 'mlapdv') || isempty(S.mlapdv)
    error('S.mlapdv is missing. Patch load_sessionPETH_all_data to store it.');
end

ML = S.mlapdv(:,1);
AP = S.mlapdv(:,2);

keep = ML >= mlRange(1) & ML < mlRange(2) & ...
       AP >= apRange(1) & AP < apRange(2);

idx = find(keep);

Sbin = S;

nTypes = numel(S.pethTypes);
for s = 1:nTypes
    Sbin.Diff_even{s} = S.Diff_even{s}(idx,:);
    Sbin.Diff_odd{s}  = S.Diff_odd{s}(idx,:);

    if isfield(S, 'MeanPos_even')
        Sbin.MeanPos_even{s} = S.MeanPos_even{s}(idx,:);
        Sbin.MeanNeg_even{s} = S.MeanNeg_even{s}(idx,:);
    end
    if isfield(S, 'MeanPos_odd')
        Sbin.MeanPos_odd{s} = S.MeanPos_odd{s}(idx,:);
        Sbin.MeanNeg_odd{s} = S.MeanNeg_odd{s}(idx,:);
    end

    Sbin.sortMetric{s} = S.sortMetric{s}(idx);
end

Sbin.globalROI  = S.globalROI(idx);
Sbin.chronicUID = S.chronicUID(idx);
Sbin.cellScore  = S.cellScore(idx);
Sbin.brainIds   = S.brainIds(idx);
Sbin.fovLabel   = S.fovLabel(idx);
Sbin.isChronic  = S.isChronic(idx);
Sbin.mlapdv     = S.mlapdv(idx,:);

if isfield(S, 'roiScale') && ~isempty(S.roiScale)
    Sbin.roiScale = S.roiScale(idx);
end
if isfield(S, 'sharedROI')
    Sbin.sharedROI = S.sharedROI(idx);
end
end