function out = mpci_tuning_heatmaps_summary(sessionPaths, alpha, opt)
% mpci_tuning_heatmaps_summary
% Collect per-session region fractions across many subjects and plot
% bias-aligned summary heatmaps (one per test).
%
% Inputs
%   sessionPaths : cellstr/string list like:
%       'Y:\Subjects\SP044\2023-08-03\002'
%   alpha        : p threshold for two-tailed criterion
%   opt.dRange   : relative day range vector (default auto from data)
%   opt.minSessionsPerBin : mask bins with fewer contributing sessions (default 3)
%   opt.weighted : true -> weighted mean by nValid, false -> simple mean (default true)
%   opt.st / opt.stPath / opt.region : passed to expand_region_ids + labeling/sorting
%   opt.biasDateFcn : function handle biasDateFcn(subjStr) -> datetime (default @getSubjectDate)
%
% Output
%   out.bigT, out.tests, out.regionList, out.regionLabels, out.dRange, out.M{ti}, out.nSess{ti}, out.nTot{ti}

if nargin < 2 || isempty(alpha), alpha = 0.05; end
if nargin < 3, opt = struct(); end

% Defaults
if ~isfield(opt,'minSessionsPerBin') || isempty(opt.minSessionsPerBin), opt.minSessionsPerBin = 3; end
if ~isfield(opt,'sumStat') || isempty(opt.sumStat), opt.sumStat = 'fracSig'; end
if ~isfield(opt,'weighted') || isempty(opt.weighted), opt.weighted = true; end
if ~isfield(opt,'region'), opt.region = {}; end
if ~isfield(opt,'st'),     opt.st = ''; end
if ~isfield(opt,'stPath'), opt.stPath = "C:\Users\Samuel\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv"; end
if ~isfield(opt,'responsiveOnly') || isempty(opt.responsiveOnly), opt.responsiveOnly = false; end
if ~isfield(opt,'responsiveAlpha') || isempty(opt.responsiveAlpha), opt.responsiveAlpha = 0.005; end
if ~isfield(opt,'biasDateFcn') || isempty(opt.biasDateFcn)
    opt.biasDateFcn = @getSubjectDate;
end
if ~isfield(opt,'ignoreGaps') || isempty(opt.ignoreGaps)
    opt.ignoreGaps = false;
end

sessionPaths = string(sessionPaths(:));
bigT = table();

% Collect all sessions into one long table
for i = 1:numel(sessionPaths)
    %Ts = mpci_tuning_FracByRegion(sessionPaths(i), alpha, opt);
    Ts = mpci_tuning_StatsByRegion(sessionPaths(i), alpha, opt);
    if ~isempty(Ts)
        bigT = [bigT; Ts]; %#ok<AGROW>
    end
end

if isempty(bigT)
    error('No data collected across sessions.');
end

% --- Optional: ignore gaps between days by using recording-order index per subject ---
% dRelIdx = (session order index) - (bias session index), per subject.
bigT.dRelIdx = nan(height(bigT),1);

subs = unique(string(bigT.subjStr), 'stable');

for si = 1:numel(subs)
    s = subs(si);
    m = string(bigT.subjStr) == s;

    days_s = bigT.day(m);
    uDays = unique(days_s);
    uDays = sort(uDays);

    % Map each session day -> ordinal index 1..K
    [~, ord] = ismember(days_s, uDays);  % ord is per-row

    % Find bias day index (may not be exactly present among uDays)
    biasDay = opt.biasDateFcn(s);
    [tfBias, biasOrd] = ismember(biasDay, uDays);

    if ~tfBias
        % If biasDay isn't exactly a recorded session day, anchor to the first day AFTER bias,
        % otherwise use last day BEFORE bias. This prevents NaNs.
        after = find(uDays >= biasDay, 1, 'first');
        if ~isempty(after)
            biasOrd = after;
        else
            biasOrd = numel(uDays);
        end
    end

    bigT.dRelIdx(m) = ord - biasOrd;
end

% Choose the x-axis variable (calendar days vs recording index)
if opt.ignoreGaps
    xVar = "dRelIdx";
else
    xVar = "dRelBias";
end

% Decide dRange (x-axis bins)
if ~isfield(opt,'dRange') || isempty(opt.dRange)
    dmin = floor(min(bigT.(xVar)));
    dmax = ceil(max(bigT.(xVar)));
    opt.dRange = dmin:dmax;
end
dRange = opt.dRange(:)';
nD = numel(dRange);

% Tests union
tests = unique(string(bigT.testName), 'stable');
nTests = numel(tests);

% Regions union
regionList = unique(double(bigT.regionId));
nR = numel(regionList);

% Load structure tree for sorting/labels
if isfield(opt,'st') && ~isempty(opt.st)
    st = opt.st;
elseif isfield(opt,'stPath') && strlength(string(opt.stPath)) > 0
    st = loadStructureTree(opt.stPath);
else
    error('Structure tree required for region labeling/sorting.');
end
rid_all = double(st.id);
acr_all = string(st.acronym);

% Sort regions by structure tree row index (idxRow)
idxRow = nan(size(regionList));
regionLabels = strings(size(regionList));
for i = 1:numel(regionList)
    ix = find(rid_all == regionList(i), 1);
    idxRow(i) = ix;
    if ~isempty(ix)
        regionLabels(i) = acr_all(ix-1);
    else
        regionLabels(i) = string(regionList(i));
    end
end
idxSort = idxRow; idxSort(isnan(idxSort)) = inf;
[~, sortIdx] = sort(idxSort);

regionList = regionList(sortIdx);
regionLabels = regionLabels(sortIdx);
idxRow = idxRow(sortIdx);

% Preallocate outputs
Mcell = cell(nTests,1);
nSessCell = cell(nTests,1);
nTotCell = cell(nTests,1);

% Build summary matrices per test
for ti = 1:nTests
    tname = tests(ti);
    Tt = bigT(string(bigT.testName) == tname, :);

    M = nan(nR, nD);
    nSess = zeros(nR, nD);
    nTot  = zeros(nR, nD);

    for ri = 1:nR
        rID = regionList(ri);
        Tr = Tt(double(Tt.regionId) == rID, :);
        if isempty(Tr), continue; end

        for di = 1:nD
            d = dRange(di);
            Td = Tr(Tr.(xVar) == d, :);
            if isempty(Td), continue; end

            x = Td.(opt.sumStat);
            w = double(Td.nValid);
            ok = ~isnan(x) & w > 0;

            if ~any(ok), continue; end

            if opt.weighted
                M(ri,di) = sum(w(ok).*x(ok)) / sum(w(ok));
            else
                M(ri,di) = mean(x(ok));
            end

            nSess(ri,di) = sum(ok);
            nTot(ri,di)  = sum(w(ok));
        end
    end

    % Mask low-support bins
    if opt.minSessionsPerBin > 1
        M(nSess < opt.minSessionsPerBin) = NaN;
    end

    Mcell{ti} = M;
    nSessCell{ti} = nSess;
    nTotCell{ti} = nTot;

    % Plot
    figure('Color','w', 'Name', char(tname), 'Position', [100, 200, 450, max(300, 15*nR + 120)]);
    
    set(gcf, 'Renderer', 'opengl');     % Ensure transparency works
    
    %imagesc(M, [0 0.5]);
    hImg = imagesc(M);
    
    if strncmp(opt.sumStat,'median',6) && strncmp(tname,'ccu',3)
        colormap(brewermap([],'*RdBu'));
        maxval = max(abs([min(M,[],'all'),max(M,[],'all')]-0.5));
        if maxval>0
            clim([-maxval,maxval]+0.5);
        end
    elseif strncmp(opt.sumStat,'median',6) && strncmp(tname,'ccMI',4)
        colormap(brewermap([],'*RdBu'));
        maxval = max(abs([min(M,[],'all'),max(M,[],'all')]));
        if maxval>0
            clim([-maxval,maxval]);
        end
    else
        colormap(parula(256));
    end
    colorbar;
    
    % Make NaNs transparent
    alphaMask = ~isnan(M);
    set(hImg, 'AlphaData', double(alphaMask));

    ax = gca;
    ax.Color = [0.85 0.85 0.85];   % light gray background
    ax.YTick = 1:nR;
    ax.YTickLabel = regionLabels;
    ax.XTick = 1:nD;
    ax.XTickLabel = string(dRange);
    ax.XTickLabelRotation = 0;

    if opt.ignoreGaps
        xlabel('Session index relative to bias (0 = bias session)');
    else
        xlabel('Days relative to bias (0 = firstBias)');
    end

    ylabel('Region');
    title(sprintf('Mean %s (%s)\n %s mean | minSessions/bin=%d', ...
        opt.sumStat, tname, ternary(opt.weighted,'weighted','unweighted'), opt.minSessionsPerBin), ...
        'Interpreter','none');

    % Mark bias day (0)
    x0 = find(dRange == 0, 1);
    if ~isempty(x0)
        xline(x0, 'r--', 'LineWidth', 2);
    end
end

out = struct();
out.bigT = bigT;
out.tests = tests;
out.regionList = regionList;
out.regionLabels = regionLabels;
out.idxRow = idxRow;
out.dRange = dRange;
out.M = Mcell;
out.nSess = nSessCell;
out.nTot = nTotCell;
out.alpha = alpha;
out.opt = opt;

end

% tiny helper
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end