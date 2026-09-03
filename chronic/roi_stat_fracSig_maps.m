function out = roi_stat_fracSig_maps(parquetFile, binSize, minROIs)
% Compute per-subject maps of fraction significant ROIs per spatial bin.
%
% A ROI is counted only if it has a non-empty cUID.
% Fraction significant in a bin = mean(h) across ROIs in that bin.
% These per-session bin fractions are then averaged across sessions
% within each subject.
%
% Inputs
% ------
% parquetFile : parquet file path
% binSize     : spatial bin size in microns, e.g. 100
% minROIs     : minimum number of ROIs per bin within a session to include
%
% Output
% ------
% out.subjects
% out.xEdges, out.yEdges
% out.meanFracSig{isub}   : mean fraction significant across sessions
% out.nDays{isub}         : number of sessions contributing to each bin
% out.meanNROIs{isub}     : mean number of ROIs per contributing session/bin

if nargin < 2 || isempty(binSize), binSize = 150; end
if nargin < 3 || isempty(minROIs), minROIs = 20; end

T = parquetread(parquetFile);

% Keep only rows with tracked IDs
hasUID = ~(ismissing(T.cUID) | T.cUID == "");
T = T(hasUID, :);

% Types
T.subject = string(T.subject);
T.date    = string(T.date);
T.session   = string(T.session);
T.cUID    = string(T.cUID);

T.ML = single(T.ML);
T.AP = single(T.AP);
T.h  = logical(T.h);

% Session key
sessionKey = T.subject + "|" + T.date + "|" + T.session;

% Global bin edges
xMin = floor(double(min(T.ML)) / binSize) * binSize;
xMax = ceil( double(max(T.ML)) / binSize) * binSize;
yMin = floor(double(min(T.AP)) / binSize) * binSize;
yMax = ceil( double(max(T.AP)) / binSize) * binSize;

xEdges = xMin:binSize:xMax;
yEdges = yMin:binSize:yMax;

nX = numel(xEdges) - 1;
nY = numel(yEdges) - 1;

subjects = unique(T.subject, 'stable');
nSub = numel(subjects);

meanFracSig = cell(nSub,1);
nDays       = cell(nSub,1);
meanNROIs   = cell(nSub,1);

for isub = 1:nSub
    subj = subjects(isub);
    fprintf('Processing subject %s (%d/%d)\n', subj, isub, nSub);

    idxSub = T.subject == subj;
    Ts = T(idxSub, :);
    Ks = sessionKey(idxSub);

    sessList = unique(Ks, 'stable');
    nSess = numel(sessList);

    sumFrac   = zeros(nY, nX, 'single');
    countFrac = zeros(nY, nX, 'uint16');
    sumN      = zeros(nY, nX, 'single');

    for n = 1:nSess
        ii = Ks == sessList(n);
        A = Ts(ii, {'ML','AP','h'});

        if isempty(A)
            continue
        end

        good = isfinite(A.ML) & isfinite(A.AP);
        if ~any(good)
            continue
        end

        ML = single(A.ML(good));
        AP = single(A.AP(good));
        h  = logical(A.h(good));

        xBin = discretize(double(ML), xEdges);
        yBin = discretize(double(AP), yEdges);

        goodBin = ~isnan(xBin) & ~isnan(yBin);
        if ~any(goodBin)
            continue
        end

        xBin = xBin(goodBin);
        yBin = yBin(goodBin);
        h    = h(goodBin);

        lin = sub2ind([nY, nX], yBin, xBin);
        [uLin, ~, g] = unique(lin);

        for ib = 1:numel(uLin)
            jj = (g == ib);
            nHere = sum(jj);

            if nHere < minROIs
                continue
            end

            fracHere = mean(single(h(jj)));

            sumFrac(uLin(ib))   = sumFrac(uLin(ib)) + fracHere;
            countFrac(uLin(ib)) = countFrac(uLin(ib)) + 1;
            sumN(uLin(ib))      = sumN(uLin(ib)) + nHere;
        end
    end

    M = nan(nY, nX, 'single');
    N = nan(nY, nX, 'single');

    valid = countFrac > 0;
    M(valid) = sumFrac(valid) ./ single(countFrac(valid));
    N(valid) = sumN(valid)    ./ single(countFrac(valid));

    meanFracSig{isub} = M;
    nDays{isub}       = countFrac;
    meanNROIs{isub}   = N;
end

out = struct();
out.subjects    = subjects;
out.xEdges      = xEdges;
out.yEdges      = yEdges;
out.meanFracSig = meanFracSig;
out.nDays       = nDays;
out.meanNROIs   = meanNROIs;
out.binSize     = binSize;
out.minROIs     = minROIs;
end