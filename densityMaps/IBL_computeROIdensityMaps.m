function D = IBL_computeROIdensityMaps(AllROIs)
%IBL_COMPUTEROIDENSITYMAPS Compute spatial maps of ROI statistics.
%
% D = IBL_computeROIdensityMaps(AllROIs)
%
% Aggregates ROI-level measurements into ML-by-AP spatial bins. The
% function computes:
%
%   - Fraction of significantly tuned ROIs
%   - Fraction of significantly positive/negative ROIs
%   - Quantile maps of the task statistic
%   - Mean and quantile maps of passive-movie cross-repeat correlation
%   - Smoothed ROI-density and session-support maps
%
% Input
% -----
% AllROIs
%   Table or structure containing, at minimum:
%
%       subject
%       date
%       session
%       p
%       h
%       stat
%
%   Spatial coordinates must be supplied either as:
%
%       ML, AP
%
%   or:
%
%       pos
%
%   where each entry in pos contains at least [ML, AP].
%
%   If present, passive-movie reliability is read from:
%
%       passiveMovieCorr
%
%   Missing passiveMovieCorr values are ignored when computing the
%   passive-movie maps.
%
% Output
% ------
% D
%   Structure containing spatial bin edges, fraction maps, task-statistic
%   quantiles, passive-movie correlation maps, density maps, plotting
%   limits, and the Allen atlas top-down projection.
%
% Notes
% -----
% Spatial maps are smoothed using a Gaussian filter. Quantile and mean maps
% are smoothed using normalized convolution so that missing bins do not
% contribute as zeros.


%% Parameters

binsize = 100;
sigma_bins = 1;

min_rois_perBin = 10;
min_rois_perBin_perSess = 5; %#ok<NASGU>
min_sess_perBin = 2;

qvec = 0.05:0.05:0.95;
nq = numel(qvec);


%% Convert input to table

if ~istable(AllROIs)
    T_full = struct2table(AllROIs);
else
    T_full = AllROIs;
end


%% Extract spatial coordinates

if ismember('pos', T_full.Properties.VariableNames)

    valid = false(height(T_full),1);

    if iscell(T_full.pos)
        valid = cellfun(@(c) ...
            isnumeric(c) && isvector(c) && numel(c) >= 2, ...
            T_full.pos);
    else
        valid = true(height(T_full),1);
    end

    T = T_full(valid,:);

    if iscell(T.pos)
        pos = cell2mat(cellfun( ...
            @(x) reshape(x(1:2),1,2), ...
            T.pos, ...
            'UniformOutput',false));
    else
        pos = T.pos(:,1:2);
    end

    ml = pos(:,1);
    ap = pos(:,2);

else
    T = T_full;

    requiredPositionVars = {'ML','AP'};
    assert( ...
        all(ismember(requiredPositionVars,T.Properties.VariableNames)), ...
        'Input must contain either pos or both ML and AP.');

    ml = T.ML;
    ap = T.AP;
end

xy = double([ml,ap]);


%% Extract ROI-level variables

requiredVars = {'subject','date','session','p','h','stat'};
missingVars = setdiff(requiredVars,T.Properties.VariableNames);

if ~isempty(missingVars)
    error( ...
        'Input is missing required variables: %s', ...
        strjoin(missingVars,', '));
end

pval = double(T.p(:));
sig = logical(T.h(:));
tstat = double(T.stat(:));

sig_hi = pval > 0.975;
sig_lo = pval < 0.025;

% Passive-movie correlation is optional
if ismember('passiveMovieCorr',T.Properties.VariableNames)
    passiveMovieCorr = double(T.passiveMovieCorr(:));
else
    passiveMovieCorr = nan(height(T),1);
end


%% Define spatial bins

validXY = all(isfinite(xy),2);

if ~any(validXY)
    error('No finite spatial coordinates found.');
end

lo = prctile(xy(validXY,:),0.1,1);
hi = prctile(xy(validXY,:),99.9,1);

xedges = lo(1):binsize:hi(1);
yedges = lo(2):binsize:hi(2);

% Ensure at least one complete bin
if numel(xedges) < 2
    xedges = [lo(1),lo(1)+binsize];
end

if numel(yedges) < 2
    yedges = [lo(2),lo(2)+binsize];
end

nbx = numel(xedges)-1;
nby = numel(yedges)-1;


%% Group ROIs by session

[G,subjlist,datelist,sesslist] = findgroups( ...
    T.subject,T.date,T.session);

ns = numel(subjlist);


%% Compute density and significance maps

SA = zeros(nbx,nby);
SS = zeros(nbx,nby);
SL = zeros(nbx,nby);
SH = zeros(nbx,nby);

n_sess_perBin = zeros(nbx,nby);

for s = 1:ns

    idxSess = G == s & validXY;

    ca = histcounts2( ...
        xy(idxSess,1), ...
        xy(idxSess,2), ...
        xedges,yedges);

    cs = histcounts2( ...
        xy(idxSess & sig,1), ...
        xy(idxSess & sig,2), ...
        xedges,yedges);

    cl = histcounts2( ...
        xy(idxSess & sig_lo,1), ...
        xy(idxSess & sig_lo,2), ...
        xedges,yedges);

    ch = histcounts2( ...
        xy(idxSess & sig_hi,1), ...
        xy(idxSess & sig_hi,2), ...
        xedges,yedges);

    n_sess_perBin = n_sess_perBin + (ca > 0);

    SA = SA + imgaussfilt( ...
        ca,sigma_bins,'FilterDomain','spatial');

    SS = SS + imgaussfilt( ...
        cs,sigma_bins,'FilterDomain','spatial');

    SL = SL + imgaussfilt( ...
        cl,sigma_bins,'FilterDomain','spatial');

    SH = SH + imgaussfilt( ...
        ch,sigma_bins,'FilterDomain','spatial');
end


%% Compute significant fractions

ratio = SS ./ SA;
ratio_lo = SL ./ SA;
ratio_hi = SH ./ SA;

insufficientROIs = SA < min_rois_perBin;
insufficientSessions = n_sess_perBin < min_sess_perBin;

invalidFractionBin = ...
    SA == 0 | insufficientROIs | insufficientSessions;

ratio(invalidFractionBin) = NaN;
ratio_lo(invalidFractionBin) = NaN;
ratio_hi(invalidFractionBin) = NaN;

mask_sess = ~insufficientSessions;


%% Assign each ROI to a spatial bin

[~,~,xb] = histcounts(xy(:,1),xedges);
[~,~,yb] = histcounts(xy(:,2),yedges);

validBin = ...
    validXY & ...
    xb > 0 & xb <= nbx & ...
    yb > 0 & yb <= nby;


%% Task-statistic quantile maps

Qstat = nan(nbx,nby,nq);

validStat = validBin & isfinite(tstat);

for ix = 1:nbx
    for iy = 1:nby

        idxBin = ...
            validStat & ...
            xb == ix & ...
            yb == iy;

        if sum(idxBin) >= min_rois_perBin
            Qstat(ix,iy,:) = quantile(tstat(idxBin),qvec);
        end
    end
end

Qstat = smoothMapStack(Qstat,sigma_bins,mask_sess);


%% Passive-movie correlation maps

QpassiveMovieCorr = nan(nbx,nby,nq);
meanPassiveMovieCorr = nan(nbx,nby);
nPassiveMovieCorr = zeros(nbx,nby);

validMovie = validBin & isfinite(passiveMovieCorr);

for ix = 1:nbx
    for iy = 1:nby

        idxBin = ...
            validMovie & ...
            xb == ix & ...
            yb == iy;

        nThisBin = sum(idxBin);
        nPassiveMovieCorr(ix,iy) = nThisBin;

        if nThisBin >= min_rois_perBin
            values = passiveMovieCorr(idxBin);

            meanPassiveMovieCorr(ix,iy) = mean(values);
            QpassiveMovieCorr(ix,iy,:) = quantile(values,qvec);
        end
    end
end

meanPassiveMovieCorr = smoothMapNormalized( ...
    meanPassiveMovieCorr,sigma_bins);

meanPassiveMovieCorr(~mask_sess) = NaN;

QpassiveMovieCorr = smoothMapStack( ...
    QpassiveMovieCorr,sigma_bins,mask_sess);


%% Output structure

D = struct();

D.xedges = xedges;
D.yedges = yedges;

% Significant fractions
D.ratio = ratio;
D.ratio_lo = ratio_lo;
D.ratio_hi = ratio_hi;

% Task statistic
D.qvec = qvec;
D.stat_quantiles = Qstat;

% Passive-movie reliability
D.passiveMovieCorr_mean = meanPassiveMovieCorr;
D.passiveMovieCorr_quantiles = QpassiveMovieCorr;
D.n_passiveMovieCorr_perBin = nPassiveMovieCorr;

% Density and support
D.SA = SA;
D.n_sess_perBin = n_sess_perBin;

% Plotting defaults
D.frac_clim = [0,0.4];
D.frac_hi_clim = [0,0.4];
D.frac_lo_clim = [0,0.4];

D.passiveMovieCorr_clim = [0,.25];

D.v_clim_sess = [1,20];
D.v_gamma = 1;

% Atlas
D.bas = aratopdown.atlas.build_topdown;

end


function maps = smoothMapStack(maps,sigmaBins,sessionMask)
%SMOOTHMAPSTACK Apply normalized Gaussian smoothing to each map layer.

nMaps = size(maps,3);

for iMap = 1:nMaps
    tmp = maps(:,:,iMap);
    tmp = smoothMapNormalized(tmp,sigmaBins);
    tmp(~sessionMask) = NaN;
    maps(:,:,iMap) = tmp;
end

end


function smoothed = smoothMapNormalized(values,sigmaBins)
%SMOOTHMAPNORMALIZED Smooth finite values without treating NaNs as zeros.

finiteMask = isfinite(values);

values0 = values;
values0(~finiteMask) = 0;

smoothedValues = imgaussfilt( ...
    values0,sigmaBins,'FilterDomain','spatial');

smoothedWeights = imgaussfilt( ...
    double(finiteMask),sigmaBins,'FilterDomain','spatial');

smoothed = smoothedValues ./ smoothedWeights;
smoothed(smoothedWeights == 0) = NaN;

end