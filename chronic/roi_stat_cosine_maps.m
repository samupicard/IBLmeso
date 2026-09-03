function out = roi_stat_cosine_maps(parquetFile, binSize, minShared, maxLag)
% ROI cosine-similarity maps using signed lags.
%
% For each subject and anchor day d, compare:
%   stat_odd(day d) vs stat_even(day d+k)
% for signed lags k = -maxLag : +maxLag.
%
% k = 0 is within-day.
%
% Required parquet columns:
%   subject, date, session, cUID, ML, AP, stat_even, stat_odd
%
% Rows without cUID are excluded.
% stat_even/stat_odd are centered by subtracting 0.5.
%
% Shuffle controls (nShuff = 200):
%   1) binShuffle*: shuffle matched ROIs within-bin
%   2) dayShuffle*: keep odd(day d) fixed, compare to even(random day)
%
% Output fields:
%   out.subjects
%   out.lags
%   out.xEdges, out.yEdges
%   out.meanCos{isub,ik}
%   out.nPairs{isub,ik}
%   out.nROIs{isub,ik}
%   out.binShuffleMean{isub,ik}
%   out.binShuffleStd{isub,ik}
%   out.binShuffleP{isub,ik}
%   out.binShuffleZ{isub,ik}
%   out.dayShuffleMean{isub,ik}
%   out.dayShuffleStd{isub,ik}
%   out.dayShuffleP{isub,ik}
%   out.dayShuffleZ{isub,ik}
%   out.nShuff

if nargin < 2 || isempty(binSize),   binSize = 150; end
if nargin < 3 || isempty(minShared), minShared = 20; end
if nargin < 4, maxLag = []; end

nShuff = 200;

T = parquetread(parquetFile);

% Keep only rows with tracked IDs
hasUID = ~(ismissing(T.cUID) | T.cUID == "");
T = T(hasUID, :);

% Types
T.subject   = string(T.subject);
T.date      = string(T.date);
T.session   = string(T.session);
T.cUID      = string(T.cUID);
T.ML        = single(T.ML);
T.AP        = single(T.AP);
T.stat_even = single(T.stat_even) - 0.5;
T.stat_odd  = single(T.stat_odd)  - 0.5;

% Parse dates
try
    dt = datetime(T.date, 'InputFormat', 'yyyy-MM-dd');
catch
    dt = datetime(T.date, 'InputFormat', 'yyyy-dd-MM');
end

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
nBins = nX * nY;

subjects = unique(T.subject, 'stable');
nSub = numel(subjects);

% Per-subject ordered session list
subjectSessions = cell(nSub,1);
for isub = 1:nSub
    idx = T.subject == subjects(isub);

    S = table(sessionKey(idx), dt(idx), T.session(idx), ...
        'VariableNames', {'sessionKey','dt','session'});

    [uSess, ia] = unique(S.sessionKey, 'stable');
    U = table(uSess, S.dt(ia), S.session(ia), ...
        'VariableNames', {'sessionKey','dt','session'});

    [~, ord] = sortrows([datenum(U.dt), str2double(U.session)]);
    subjectSessions{isub} = U.sessionKey(ord);
end

maxPossibleLag = max(cellfun(@numel, subjectSessions)) - 1;
if isempty(maxLag)
    maxLagUsed = maxPossibleLag;
else
    maxLagUsed = min(maxLag, maxPossibleLag);
end

lags = -maxLagUsed:maxLagUsed;
nLag = numel(lags);

meanCos = cell(nSub, nLag);
nPairs  = cell(nSub, nLag);
nROIs   = cell(nSub, nLag);

binShuffleCos = cell(nSub, nLag);
binShuffleP   = cell(nSub, nLag);

dayShuffleCos = cell(nSub, nLag);
dayShuffleP   = cell(nSub, nLag);

for isub = 1:nSub
    subj = subjects(isub);
    fprintf('Processing subject %s (%d/%d)\n', subj, isub, nSub);

    idxSub = T.subject == subj;
    Ts = T(idxSub, :);
    Ks = sessionKey(idxSub);

    sessList = subjectSessions{isub};
    nSess = numel(sessList);

    % Store each session as plain arrays, sorted by cUID
    sess = repmat(struct( ...
        'cUID', strings(0,1), ...
        'ML', single([]), ...
        'AP', single([]), ...
        'odd', single([]), ...
        'even', single([])), nSess, 1);

    for i = 1:nSess
        ii = Ks == sessList(i);
        X = Ts(ii, {'cUID','ML','AP','stat_odd','stat_even'});

        % One row per cUID per session
        [uID, ia] = unique(X.cUID, 'stable');
        X = X(ia, :);
        X.cUID = uID;

        % Sort by cUID for faster intersect behavior
        [cSorted, ord] = sort(X.cUID);
        sess(i).cUID = cSorted;
        sess(i).ML   = single(X.ML(ord));
        sess(i).AP   = single(X.AP(ord));
        sess(i).odd  = single(X.stat_odd(ord));
        sess(i).even = single(X.stat_even(ord));
    end

    % Accumulators per lag
    sumCos   = cell(nLag,1);
    countCos = cell(nLag,1);
    sumN     = cell(nLag,1);

    sumBinP   = cell(nLag,1);
    sumDayP   = cell(nLag,1);

    sumBinCos = cell(nLag,1);   % each is nBins x nShuff
    cntBinCos = cell(nLag,1);   % each is nBins x nShuff

    sumDayCos = cell(nLag,1);   % each is nBins x nShuff
    cntDayCos = cell(nLag,1);   % each is nBins x nShuff

    for ik = 1:nLag
        sumCos{ik}   = zeros(nBins,1,'single');
        countCos{ik} = zeros(nBins,1,'uint16');
        sumN{ik}     = zeros(nBins,1,'single');

        sumBinP{ik}   = zeros(nBins,1,'single');
        sumDayP{ik}   = zeros(nBins,1,'single');

        sumBinCos{ik} = zeros(nBins,nShuff,'single');
        cntBinCos{ik} = zeros(nBins,nShuff,'uint16');

        sumDayCos{ik} = zeros(nBins,nShuff,'single');
        cntDayCos{ik} = zeros(nBins,nShuff,'uint16');
    end

    % Loop over anchor day d
    for d = 1:nSess
        fprintf('  day %d / %d\n', d, nSess);

        A = sess(d);
        if isempty(A.cUID)
            continue
        end

        % Cache all pairwise data for this anchor day d
        pairCache = cell(nSess,1);

        d2min = max(1, d - maxLagUsed);
        d2max = min(nSess, d + maxLagUsed);

        for d2 = d2min:d2max
            B = sess(d2);
            if isempty(B.cUID)
                continue
            end

            [~, ia, ib] = intersect(A.cUID, B.cUID, 'stable');
            if isempty(ia)
                continue
            end

            ML1 = A.ML(ia);
            AP1 = A.AP(ia);
            odd = A.odd(ia);

            ML2 = B.ML(ib);
            AP2 = B.AP(ib);
            even = B.even(ib);

            good = isfinite(ML1) & isfinite(AP1) & ...
                   isfinite(ML2) & isfinite(AP2) & ...
                   isfinite(odd) & isfinite(even);
            if ~any(good)
                continue
            end

            ML1 = ML1(good);
            AP1 = AP1(good);
            ML2 = ML2(good);
            AP2 = AP2(good);
            odd = odd(good);
            even = even(good);

            MLmid = 0.5 * (ML1 + ML2);
            APmid = 0.5 * (AP1 + AP2);

            xBin = discretize(double(MLmid), xEdges);
            yBin = discretize(double(APmid), yEdges);

            goodBin = ~isnan(xBin) & ~isnan(yBin);
            if ~any(goodBin)
                continue
            end

            xBin = xBin(goodBin);
            yBin = yBin(goodBin);
            odd  = odd(goodBin);
            even = even(goodBin);

            lin = sub2ind([nY, nX], yBin, xBin);

            pairCache{d2} = struct( ...
                'lin',  lin(:), ...
                'odd',  odd(:), ...
                'even', even(:));
        end

        % Real and bin-shuffle computations for signed lags
        for ik = 1:nLag
            k = lags(ik);
            d2 = d + k;
            if d2 < 1 || d2 > nSess
                continue
            end

            P = pairCache{d2};
            if isempty(P)
                continue
            end

            lin  = P.lin;
            odd  = P.odd;
            even = P.even;

            % Real cosine via accumarray
            cnt   = accumarray(lin, 1,         [nBins 1], @sum, 0);
            sumXY = accumarray(lin, single(odd.*even), [nBins 1], @sum, single(0));
            sumX2 = accumarray(lin, single(odd.^2),    [nBins 1], @sum, single(0));
            sumY2 = accumarray(lin, single(even.^2),   [nBins 1], @sum, single(0));

            denom = sqrt(sumX2 .* sumY2);
            cVec = sumXY ./ denom;
            cVec(cnt < minShared | denom <= 0 | ~isfinite(denom)) = NaN;

            validBins = find(isfinite(cVec));
            if isempty(validBins)
                continue
            end

            sumCos{ik}(validBins)   = sumCos{ik}(validBins)   + single(cVec(validBins));
            countCos{ik}(validBins) = countCos{ik}(validBins) + 1;
            sumN{ik}(validBins)     = sumN{ik}(validBins)     + single(cnt(validBins));

            % Group bin members once for shuffle controls
            [uLin, ~, g] = unique(lin);

            % Bin-shuffle null
            for ibin = 1:numel(uLin)
                binId = uLin(ibin);
                jj = (g == ibin);

                vOdd  = odd(jj);
                vEven = even(jj);

                [cReal, nHere] = local_bin_cosine(vOdd, vEven, minShared, false);
                if ~isfinite(cReal)
                    continue
                end

                cShuff = nan(nShuff,1,'single');
                for ish = 1:nShuff
                    cShuff(ish) = local_bin_cosine(vOdd, vEven, minShared, true);
                end
                cShuff = cShuff(isfinite(cShuff));
                if isempty(cShuff)
                    continue
                end

                p = (sum(cShuff >= cReal) + 1) / (numel(cShuff) + 1);

                if isfinite(p)
                    sumBinP{ik}(binId) = sumBinP{ik}(binId) + p;
                end

                validSh = isfinite(cShuff);
                if any(validSh)
                    shIdx = find(validSh);
                    sumBinCos{ik}(binId, shIdx) = sumBinCos{ik}(binId, shIdx) + reshape(single(cShuff(validSh)), 1, []);
                    cntBinCos{ik}(binId, shIdx) = cntBinCos{ik}(binId, shIdx) + 1;
                end
            end

            % Day-shuffle null: same odd(day d), random even(day r)
            randDays = randi(nSess, nShuff, 1);

            % Precompute all sampled random-day pair caches only once
            randCache = cell(nShuff,1);
            for ish = 1:nShuff
                randCache{ish} = pairCache{randDays(ish)};
                if isempty(randCache{ish})
                    % pairCache only covers d2 within lag window. Build on demand.
                    r = randDays(ish);
                    B = sess(r);
                    if ~isempty(B.cUID)
                        [~, ia, ib] = intersect(A.cUID, B.cUID, 'stable');
                        if ~isempty(ia)
                            ML1 = A.ML(ia);
                            AP1 = A.AP(ia);
                            oddr = A.odd(ia);
                            ML2 = B.ML(ib);
                            AP2 = B.AP(ib);
                            evenr = B.even(ib);

                            good = isfinite(ML1) & isfinite(AP1) & ...
                                   isfinite(ML2) & isfinite(AP2) & ...
                                   isfinite(oddr) & isfinite(evenr);
                            if any(good)
                                ML1 = ML1(good);
                                AP1 = AP1(good);
                                ML2 = ML2(good);
                                AP2 = AP2(good);
                                oddr = oddr(good);
                                evenr = evenr(good);

                                MLmid = 0.5 * (ML1 + ML2);
                                APmid = 0.5 * (AP1 + AP2);

                                xBin = discretize(double(MLmid), xEdges);
                                yBin = discretize(double(APmid), yEdges);

                                goodBin = ~isnan(xBin) & ~isnan(yBin);
                                if any(goodBin)
                                    xBin = xBin(goodBin);
                                    yBin = yBin(goodBin);
                                    oddr = oddr(goodBin);
                                    evenr = evenr(goodBin);

                                    randCache{ish} = struct( ...
                                        'lin',  sub2ind([nY, nX], yBin, xBin), ...
                                        'odd',  oddr(:), ...
                                        'even', evenr(:));
                                end
                            end
                        end
                    end
                end
            end

            % Use the same valid bins as the real comparison
            for ibin = 1:numel(validBins)
                binId = validBins(ibin);

                % real bin value from actual lag-k comparison
                jjReal = (lin == binId);
                [cReal, ~] = local_bin_cosine(odd(jjReal), even(jjReal), minShared, false);
                if ~isfinite(cReal)
                    continue
                end

                cShuffDay = nan(nShuff,1,'single');

                for ish = 1:nShuff
                    R = randCache{ish};
                    if isempty(R)
                        continue
                    end

                    jj = (R.lin == binId);
                    if ~any(jj)
                        continue
                    end

                    cShuffDay(ish) = local_bin_cosine(R.odd(jj), R.even(jj), minShared, false);
                end

                cShuffDay = cShuffDay(isfinite(cShuffDay));
                if isempty(cShuffDay)
                    continue
                end

                p = (sum(cShuffDay >= cReal) + 1) / (numel(cShuffDay) + 1);

                if isfinite(p)
                    sumDayP{ik}(binId) = sumDayP{ik}(binId) + p;
                end

                validSh = isfinite(cShuffDay);
                if any(validSh)
                    shIdx = find(validSh);
                    sumDayCos{ik}(binId, shIdx) = sumDayCos{ik}(binId, shIdx) + reshape(single(cShuffDay(validSh)), 1, []);
                    cntDayCos{ik}(binId, shIdx) = cntDayCos{ik}(binId, shIdx) + 1;
                end
            end
        end
    end

    % Finalize subject outputs
    for ik = 1:nLag
        valid = countCos{ik} > 0;

        M    = nan(nBins,1,'single');
        N    = nan(nBins,1,'single');
        BinP = nan(nBins,1,'single');
        DayP = nan(nBins,1,'single');

        M(valid)    = sumCos{ik}(valid) ./ single(countCos{ik}(valid));
        N(valid)    = sumN{ik}(valid)   ./ single(countCos{ik}(valid));
        BinP(valid) = sumBinP{ik}(valid) ./ single(countCos{ik}(valid));
        DayP(valid) = sumDayP{ik}(valid) ./ single(countCos{ik}(valid));

        % mean shuffle cosine for each shuffle index, per bin
        BinCos = nan(nBins, nShuff, 'single');
        DayCos = nan(nBins, nShuff, 'single');

        validBinSh = cntBinCos{ik} > 0;
        BinCos(validBinSh) = sumBinCos{ik}(validBinSh) ./ single(cntBinCos{ik}(validBinSh));

        validDaySh = cntDayCos{ik} > 0;
        DayCos(validDaySh) = sumDayCos{ik}(validDaySh) ./ single(cntDayCos{ik}(validDaySh));

        meanCos{isub,ik}      = reshape(M,    [nY nX]);
        nPairs{isub,ik}       = reshape(countCos{ik}, [nY nX]);
        nROIs{isub,ik}        = reshape(N,    [nY nX]);

        binShuffleP{isub,ik}  = reshape(BinP, [nY nX]);
        dayShuffleP{isub,ik}  = reshape(DayP, [nY nX]);

        binShuffleCos{isub,ik} = reshape(BinCos, [nY nX nShuff]);
        dayShuffleCos{isub,ik} = reshape(DayCos, [nY nX nShuff]);
    end
end

out = struct();
out.subjects = subjects;
out.lags     = lags;
out.xEdges   = xEdges;
out.yEdges   = yEdges;

out.meanCos  = meanCos;
out.nPairs   = nPairs;
out.nROIs    = nROIs;

out.binShuffleCos = binShuffleCos;
out.binShuffleP   = binShuffleP;

out.dayShuffleCos = dayShuffleCos;
out.dayShuffleP   = dayShuffleP;

out.nShuff = nShuff;
end

function [c, nHere] = local_bin_cosine(v1, v2, minShared, doShuffle)
if nargin < 4 || isempty(doShuffle)
    doShuffle = false;
end

v1 = single(v1(:));
v2 = single(v2(:));

valid = isfinite(v1) & isfinite(v2);
v1 = v1(valid);
v2 = v2(valid);

nHere = numel(v1);
if nHere < minShared
    c = NaN;
    return
end

if doShuffle
    v2 = v2(randperm(nHere));
end

denom = norm(v1) * norm(v2);
if denom <= 0 || ~isfinite(denom)
    c = NaN;
    return
end

c = dot(v1, v2) / denom;
if ~isfinite(c)
    c = NaN;
end
end