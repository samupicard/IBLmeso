function out = plot_sessionPETH(datpath, pethType, varargin)
% plot_sessionPETH
%
% Loads PETH subtype (e.g. 'stimOn_contrastDiff') for all FOVs in a session,
% filters ROIs by cellClassifier (>0.5) and optionally brain region,
% sorts neurons by odd-trial diff between extreme conditions in twin_ev,
% and plots a 2 x nConditions grid of imagesc (odd top, even bottom).
%
% Requirements:
%   readNPY (and Allen structure tree loader if using region filtering)
%
% Expected files per FOV folder:
%   mpciROIs.cellClassifier.npy            [nROIs]
%   mpciROIs.brainLocationIds.npy          [nROIs]
%   mpciROIs.PETHavgNorm_<evnt>_<field>_odd.npy   [nROIs x nConds x nTime]
%   mpciROIs.PETHavgNorm_<evnt>_<field>_even.npy  [nROIs x nConds x nTime]
%   PETHavgNorm_<evnt>_<field>.timeValues.npy [nTime]
%   PETHavgNorm_<evnt>_<field>.conditionValues.npy [nConds]
%
% Fall-back (if odd/even not present):
%   mpciROIs.PETHavgNorm_<evnt>_<field>.npy
%
% Outputs:
%   out struct with concatenated data, masks, sorting index, etc.

%% parse inputs
p = inputParser;
p.addParameter('classifierThresh', 0.5, @(x)isnumeric(x) && isscalar(x));
p.addParameter('onlyResponsive', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('onlyTuned', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('region', {}, @(x) iscell(x) && all(cellfun(@ischar, x)));
p.addParameter('stPath', 'C:\Users\Samuel\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv', @(s)ischar(s) || isstring(s));
p.addParameter('twin_ev', [], @(x)isnumeric(x) && (isempty(x) || (numel(x)==2)));
p.addParameter('T', [], @(x)isnumeric(x) && isvector(x));  % optional explicit timebase
p.addParameter('condVals', [], @(x)isnumeric(x) || iscell(x)); % optional override
p.addParameter('sortAbs', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('caxis', [], @(x)isnumeric(x) && (isempty(x) || numel(x)==2));
p.addParameter('showColorbar', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('subsampleN', [], @(x)isnumeric(x) && (isempty(x) || (isscalar(x) && x>0)));
p.addParameter('subsampleSeed', 0, @(x)isnumeric(x) && isscalar(x)); % 0 => no rng shuffle, deterministic

p.parse(varargin{:});
opt = p.Results;
opt.stPath = char(opt.stPath);
opt.sortAbs = logical(opt.sortAbs);
opt.showColorbar = logical(opt.showColorbar);
if ~isempty(opt.subsampleN)
    opt.subsampleN = round(opt.subsampleN);
end

%% parse pethType => evnt + trialTypeField
% expected: 'stimOn_contrastDiff'
parts = split(string(pethType), "_");
if numel(parts) < 2
    error('pethType must look like "stimOn_contrastDiff"');
end
evnt = char(parts(1));
trialTypeField = char(strjoin(parts(2:end), "_")); % in case field contains underscores later

%% defaults per known PETH types
[defaultTwinEv, defaultCondVals, defaultT] = get_twins(evnt, trialTypeField);

if isempty(opt.twin_ev), opt.twin_ev = defaultTwinEv; end
if isempty(opt.condVals), opt.condVals = defaultCondVals; end

% timebase: prefer event-specific T, then defaultT; can be overwritten with user-prodived T; otherwise infer indices only
if isempty(opt.T)
    T = defaultT; % may still be empty
else
    T = opt.T;
end

%% find FOV folders
alfDir = fullfile(datpath, 'alf');
if ~exist(alfDir, 'dir')
    error('Cannot find alf folder: %s', alfDir);
end

d = dir(fullfile(alfDir, 'FOV_*'));
d = d([d.isdir]);
if isempty(d)
    error('No FOV_* folders found in %s', alfDir);
end

%% load + concat across FOVs
Podd_all       = [];
Peven_all      = [];
cellScore_all  = [];
brainId_all    = [];
fovLabel_all   = [];
taskResp_all   = [];
taskTuned_all  = [];

warnedFallback = false;

for i = 1:numel(d)
    fovFolder = fullfile(d(i).folder, d(i).name);

    cellFile = fullfile(fovFolder, 'mpciROIs.cellClassifier.npy');
    brnFile  = fullfile(fovFolder, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy');
    typeFile = fullfile(fovFolder, 'mpciROIs.mpciROITypes.npy');
    respFile  = fullfile(fovFolder, 'mpciROIs.taskResponsiveP.tsv');
    tunedFile = fullfile(fovFolder, 'mpciROIs.taskTunedP.tsv');

    if ~isfile(cellFile) || ~isfile(brnFile) || ~isfile(typeFile)
        warning('Missing classifier/brainLocation/roiType files in %s (skipping FOV)', fovFolder);
        continue
    end

    cellScore = readNPY(cellFile);  cellScore = cellScore(:);
    brainIds  = readNPY(brnFile);   brainIds  = brainIds(:);
    roiTypes  = readNPY(typeFile);  roiTypes  = roiTypes(:);

    if numel(cellScore) ~= numel(brainIds) || numel(cellScore) ~= numel(roiTypes)
        warning('ROI metadata size mismatch in %s (skipping FOV)', fovFolder);
        continue
    end

    % Keep only somatic ROIs, following newer convention
    keepType = (double(roiTypes) == 1);
    if ~any(keepType)
        continue
    end

    nRois = numel(cellScore);

    % responsive filter from TSV
    taskResp = true(nRois,1);
    respNames = string.empty(1,0);
    Resp = [];

    if opt.onlyResponsive
        if isfile(respFile)
            try
                respTab = readtable(respFile, ...
                    'FileType','text', ...
                    'Delimiter','\t', ...
                    'VariableNamingRule','preserve');
                respNames = string(respTab.Properties.VariableNames);
                Resp = table2array(respTab);

                if size(Resp,1) ~= nRois
                    warning('Responsive table size mismatch in %s; using all ROIs for responsive filter.', fovFolder);
                    Resp = [];
                    respNames = string.empty(1,0);
                else
                    % Same convention as old code: any significant column after Bonferroni correction
                    taskResp = any(Resp < (0.05 / size(Resp,2)), 2);
                end
            catch ME
                warning('Could not load %s (%s). Using all ROIs for responsive filter.', respFile, ME.message);
            end
        else
            warning('Could not find mpciROIs.taskResponsiveP.tsv in %s; using all ROIs.', fovFolder);
        end
    end

    % tuned filter from TSV
    taskTuned = true(nRois,1);

    if opt.onlyTuned
        if isfile(tunedFile)
            try
                tunedTab = readtable(tunedFile, ...
                    'FileType','text', ...
                    'Delimiter','\t', ...
                    'VariableNamingRule','preserve');
                tunedNames = string(tunedTab.Properties.VariableNames);
                tunedPvals = table2array(tunedTab);

                if size(tunedPvals,1) ~= nRois
                    warning('Tuned table size mismatch in %s; using all ROIs for tuning filter.', fovFolder);
                else
                    if islogical(opt.onlyTuned)
                        % any tuned test significant, Bonferroni across columns
                        taskTuned = any(tunedPvals < (0.05 / size(tunedPvals,2)), 2);
                    else
                        % opt.onlyTuned may be numeric index, char, string, or cellstr
                        selIdx = [];

                        if isnumeric(opt.onlyTuned)
                            selIdx = opt.onlyTuned(:)';
                        elseif ischar(opt.onlyTuned) || isstring(opt.onlyTuned)
                            selNames = string(opt.onlyTuned);
                            selIdx = find(ismember(tunedNames, selNames));
                        elseif iscellstr(opt.onlyTuned) || iscell(opt.onlyTuned)
                            selNames = string(opt.onlyTuned);
                            selIdx = find(ismember(tunedNames, selNames));
                        end

                        selIdx = selIdx(selIdx >= 1 & selIdx <= size(tunedPvals,2));

                        if isempty(selIdx)
                            warning('Requested tuning selection not found in %s; using all ROIs.', fovFolder);
                        else
                            % preserve old percentile-style logic, but vectorized correctly
                            selP = tunedPvals(:, selIdx);
                            taskTuned = any((selP < 0.025) | (selP > 0.975), 2);
                        end
                    end
                end
            catch ME
                warning('Could not load %s (%s). Using all ROIs for tuning filter.', tunedFile, ME.message);
            end
        else
            warning('Could not find mpciROIs.taskTunedP.tsv in %s; using all ROIs.', fovFolder);
        end
    end

    % load PETHs
    oddFile  = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s_odd.npy', evnt, trialTypeField));
    evenFile = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s_even.npy', evnt, trialTypeField));
    avgFile  = fullfile(fovFolder, sprintf('mpciROIs.PETHavgNorm_%s_%s.npy', evnt, trialTypeField));
    timeFile = fullfile(fovFolder, sprintf('PETHavgNorm_%s_%s.timeValues.npy', evnt, trialTypeField));

    if isfile(oddFile) && isfile(evenFile)
        Podd  = readNPY(oddFile);
        Peven = readNPY(evenFile);
    elseif isfile(avgFile)
        Podd  = readNPY(avgFile);
        Peven = Podd;
        if ~warnedFallback
            warning('Odd/Even files not found. Falling back to avgNorm for both rows.');
            warnedFallback = true;
        end
    else
        warning('Missing PETH files for %s in %s (skipping FOV)', pethType, fovFolder);
        continue
    end

    if isfile(timeFile)
        T = readNPY(timeFile);
    elseif ~isempty(opt.T),
        T = opt.T;
    else
        warning('Missing time vector files for %s in %s, and no manual time vector provided (skipping FOV)', pethType, fovFolder);
        continue
    end

    % sanity checks
    if size(Podd,1) ~= nRois || size(Peven,1) ~= nRois
        warning('ROI count mismatch in %s (skipping FOV)', fovFolder);
        continue
    end

    if numel(T) ~= size(Podd,3)
        warning('Time vector mismatch in %s (skipping FOV)', fovFolder);
        continue
    end

    % Apply ROI-type filter after all size checks
    keep = keepType(:);

    Podd      = Podd(keep,:,:);
    Peven     = Peven(keep,:,:);
    cellScore = cellScore(keep);
    brainIds  = brainIds(keep);
    taskResp  = taskResp(keep);
    taskTuned = taskTuned(keep);

    % concat
    Podd_all       = cat(1, Podd_all, Podd);
    Peven_all      = cat(1, Peven_all, Peven);
    cellScore_all  = cat(1, cellScore_all, cellScore);
    brainId_all    = cat(1, brainId_all, brainIds);
    taskResp_all   = cat(1, taskResp_all, taskResp);
    taskTuned_all  = cat(1, taskTuned_all, taskTuned);

    % label fov number if parseable
    fovNum = nan;
    tok = regexp(d(i).name, 'FOV_(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        fovNum = str2double(tok{1});
    end
    fovLabel_all = cat(1, fovLabel_all, repmat(fovNum, sum(keep), 1));
end

if isempty(Podd_all)
    error('No data loaded for %s from %s', pethType, datpath);
end

nROIs  = size(Podd_all,1);
nConds = size(Podd_all,2);
nTime  = size(Podd_all,3);

%% build time vector if still empty
if isempty(T)
    % fall back: use sample index as "time"
    T = 1:nTime;
end
T = T(:)';

%% cell classifier & task-responsive filter
isCell = cellScore_all > opt.classifierThresh;
isResp = logical(taskResp_all);

%% optional region filter
isRegion = true(nROIs,1);
regionIdsExpanded = [];

if ~isempty(opt.region)

    st = loadStructureTree(opt.stPath);

    % normalize fields
    acr_all = string(st.acronym);
    rid_all = double(st.id);

    % depth column name safety
    if ismember('depth', st.Properties.VariableNames)
        depth_all = double(st.depth);
    else
        error('Structure tree table must contain a "depth" column.');
    end

    % 1) expand region "tokens" (e.g. HVAs -> {'VISa','VISl',...})
    %    and keep wildcards/acronyms as-is
    regList = {};  % cellstr of items to match (may include wildcards like 'VIS*')
    for r = 1:numel(opt.region)
        tok = opt.region{r};

        % user-defined mapping (return {} if not a token)
        acrList = regionTokenMap(tok);

        if ~isempty(acrList)
            % expand token into explicit acronyms (these are plain acronyms, no wildcards)
            regList = [regList, acrList(:)']; %#ok<AGROW>
        else
            regList = [regList, {tok}]; %#ok<AGROW>
        end
    end

    % 2) build the union of matching IDs
    regionIdsExpanded = [];

    for r = 1:numel(regList)
        reg = char(regList{r});

        % wildcard case: e.g. 'VIS*' -> any acronym starting with 'VIS'
        if endsWith(reg, '*')
            prefix = reg(1:end-1);
            match = startsWith(acr_all, string(prefix));
            regionIdsExpanded = [regionIdsExpanded; rid_all(match)]; %#ok<AGROW>
            continue
        end

        % plain acronym case: 'VISa'
        % exact match always included
        exact = (acr_all == string(reg));
        if ~any(exact)
            % nothing matched; skip
            continue
        end

        % include the exact region id
        regionIdsExpanded = [regionIdsExpanded; rid_all(exact)]; %#ok<AGROW>

        % depth-gated prefix matching:
        % only if the matched acronym has depth >= 7
        d = depth_all(find(exact, 1, 'first'));
        if d >= 7
            match = startsWith(acr_all, string(reg));
            regionIdsExpanded = [regionIdsExpanded; rid_all(match)]; %#ok<AGROW>
        end
    end

    regionIdsExpanded = unique(regionIdsExpanded);

    if isempty(regionIdsExpanded)
        warning('None of the requested regions/tokens matched: %s', strjoin(opt.region, ', '));
        isRegion = false(nROIs,1);
    else
        isRegion = ismember(double(brainId_all), regionIdsExpanded);
    end
    
    if numel(opt.region) == 1
        regs = opt.region{1};
    else
        regs = strjoin(opt.region, ', ');
    end
else
    regs = {'All'};
end

keep = isCell & isRegion & isResp;
if ~any(keep)
    error('No ROIs left after filtering (cellThresh=%.2f, region="%s, taskResp=")', opt.classifierThresh, opt.region, sum(isResp));
end

Podd_keep  = Podd_all(keep,:,:);
Peven_keep = Peven_all(keep,:,:);

%% sorting index based on odd-trial response-difference in twin_ev window
twin = opt.twin_ev;
if isempty(twin) || numel(twin)~=2
    error('twin_ev must be a 1x2 vector like [0 0.4]');
end

% indices of time window
tMask = (T >= twin(1)) & (T <= twin(2));
if ~any(tMask)
    % if T is just indices, interpret twin as indices (best effort)
    if isnumeric(T) && isequal(T, 1:nTime)
        tMask = (T >= twin(1)) & (T <= twin(2));
    end
end
if ~any(tMask)
    error('twin_ev [%.3f %.3f] does not overlap timebase.', twin(1), twin(2));
end

respLow  = squeeze(mean(Podd_keep(:,1,tMask), 3));      % nNeurons x 1
respHigh = squeeze(mean(Podd_keep(:,end,tMask), 3));    % nNeurons x 1
dResp = respHigh - respLow;
miResp = (respHigh - respLow) / (respHigh + respLow);

if opt.sortAbs
    [~, sortIdx] = sort(abs(dResp), 'descend');
else
    [~, sortIdx] = sort(dResp, 'descend');
end

Podd_sorted  = Podd_keep(sortIdx,:,:);
Peven_sorted = Peven_keep(sortIdx,:,:);
Deven_sorted = squeeze(Peven_keep(sortIdx,end,:) - Peven_keep(sortIdx,1,:));

% optional subsample neurons so rows are visible
if ~isempty(opt.subsampleN) && opt.subsampleN < size(Podd_sorted,1)
    nAvail = size(Podd_sorted,1);

    if opt.subsampleSeed ~= 0
        rng(opt.subsampleSeed);
        sel = sort(randperm(nAvail, opt.subsampleN)); % keep ascending so order is stable-ish
    else
        % deterministic: take evenly spaced indices through sorted list
        sel = unique(round(linspace(1, nAvail, opt.subsampleN)));
        % if rounding reduced count, pad by adding neighbors
        while numel(sel) < opt.subsampleN
            cand = setdiff(1:nAvail, sel);
            sel = sort([sel, cand(1)]);
        end
    end

    Podd_sorted  = Podd_sorted(sel,:,:);
    Peven_sorted = Peven_sorted(sel,:,:);

    % keep track of which sorted neurons survived
    sortIdx = sortIdx(sel);
end


%% plot
figure('Position',[1000,50,100+nConds*50,920],'Color','w',...
    'Name',sprintf('%s | %s', sessNameFromPath(datpath), pethType));
%tiledlayout(2, nConds, 'TileSpacing', 'compact', 'Padding', 'compact');
tiledlayout(1, nConds+1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Determine condition labels
condVals = opt.condVals;
if isnumeric(condVals) && numel(condVals)==nConds
    condLabels = arrayfun(@(x) sprintf('%g', x), condVals, 'UniformOutput', false);
else
    % fallback: 1..nConds
    condLabels = arrayfun(@(c) sprintf('cond %d', c), 1:nConds, 'UniformOutput', false);
end

% set a single caxis across all subplots (unless user provided one)
if isempty(opt.caxis)
    % compute from all displayed data (odd+even)
    allVals = cat(1, Podd_sorted(:), Peven_sorted(:));
    % robust limits (avoid single hot pixel dominating)
    clim = prctile(allVals, [1 99.9]);
    if clim(1) == clim(2)
        clim = [min(allVals) max(allVals)];
        if clim(1) == clim(2)
            clim = [clim(1)-eps clim(2)+eps];
        end
    end
    opt.caxis = clim;
end

ax = [];
for c = 1:nConds
    % % --- top row (odd)
    % nexttile(c);
    % imgPodd = squeeze(Podd_sorted(:,c,:)); % nNeurons x nTime
    % imagesc(T, 1:size(imgPodd,1), imgPodd);
    % colormap(brewermap([],'Greys'));
    % axis tight
    % set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'TickLength', [0 0]);
    % title(sprintf('%s', condLabels{c}), 'Interpreter','none');
    % hold on; xline(0,'k-'); hold off
    % caxis(opt.caxis);    
    
    % --- bottom row (even)
    %nexttile(nConds + c);
    ax(c) = nexttile(c);
    imgPeven = squeeze(Peven_sorted(:,c,:));
    imagesc(T, 1:size(imgPeven,1), imgPeven);
    colormap(ax(c),brewermap([],'Greys'));
    axis tight
    set(ax(c), 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'TickLength', [0 0]);
    title(sprintf('%s', condLabels{c}), 'Interpreter','none');
    hold on; xline(0,'k-'); hold off
    caxis(ax(c),opt.caxis);    

end

% last column: difference plot (even)
ax2 = nexttile;
imgDeven = squeeze(Peven_sorted(:,end,:) - Peven_sorted(:,1,:)); 
imagesc(T, 1:size(imgDeven,1), imgDeven);
colormap(ax2,brewermap([],'RdBu'));
axis tight
set(ax2, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'TickLength', [0 0]);
title('DIFF', 'Interpreter','none');
hold on; xline(0,'k-'); hold off
caxis(ax2,0.8*opt.caxis(end)*[-1,1]);

if opt.showColorbar
    cb = colorbar;
    cb.Layout.Tile = 'east';
end

% ---- figure-level scalebar outside bottom-left axes ----
% (drawn as annotations in normalized figure coordinates)

% get bottom-left axes handle (row 2, col 1)
%axBL = nexttile(nConds + 1);
axBL = nexttile(1);
pos = axBL.Position;  % [x y w h] in normalized figure units

% lengths in data units
dx_sec = 1; 
dy_nrn = 10 * round(size(Peven_sorted,1)/500); 

% convert data lengths -> normalized figure lengths relative to that axes
Tspan = T(end) - T(1);
if Tspan <= 0
    Tspan = 1; % safety
end
nNeur = size(Peven_sorted,1);
if nNeur <= 0
    nNeur = 1; % safety
end

dx_norm = (dx_sec / Tspan) * pos(3);
dy_norm = (dy_nrn / nNeur) * pos(4);

% anchor point just OUTSIDE bottom-left axes:
% a little left of axes and a little below axes
marginX = 0.010;   % tweak if you want
marginY = 0.025;

x0 = pos(1) - marginX;
y0 = pos(2) - marginY;

% keep within figure bounds (avoid negative coordinates)
x0 = max(x0, 0.02);
y0 = max(y0, 0.02);

% draw horizontal bar (1s)
annotation(gcf, 'line', [x0, x0 + dx_norm], [y0, y0], ...
    'Color', 'k', 'LineWidth', 2);

% draw vertical bar (N neurons) - upwards from y0
annotation(gcf, 'line', [x0, x0], [y0, y0 + dy_norm], ...
    'Color', 'k', 'LineWidth', 2);

% labels
% annotation(gcf, 'textbox', [x0 + dx_norm/2 - 0.03, y0 - 0.03, 0.06, 0.03], ...
%     'String', sprintf('%ds',dx_sec), 'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
%     'Color', 'k', 'FontSize', 8);
% 
% annotation(gcf, 'textbox', [x0 - 0.03, y0 + dy_norm/2 - 0.015, 0.05, 0.03], ...
%     'String', sprintf('%d',dy_nrn), 'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
%     'Color', 'k', 'FontSize', 8);


% add grand title
sgtitle(sprintf('%s-aligned\nsplit by %s\n%s (n=%d)', ...
    evnt, trialTypeField, regs, size(Podd_sorted,1)),...
    'Interpreter','none');
%sgtitle(sprintf('%s | sort by DIFF in twin=[%.2f %.2f] | n=%d', ...
%    pethType, twin(1), twin(2), size(Podd_sorted,1)), 'Interpreter','none');

%% outputs
out = struct();
out.datpath = datpath;
out.pethType = pethType;
out.evnt = evnt;
out.trialTypeField = trialTypeField;
out.T = T;
out.twin_ev = twin;

out.nLoadedROIs = nROIs;
out.keepMask = keep;
out.sortIdx = sortIdx;

out.cellScore = cellScore_all;
out.brainIds = brainId_all;
out.fovLabel = fovLabel_all;
out.region = opt.region;
out.regionIdsExpanded = regionIdsExpanded;

out.Podd_sorted = Podd_sorted;
out.Peven_sorted = Peven_sorted;
out.MI = miResp;
out.diff = dResp;

end

function nm = sessNameFromPath(datpath)
% best-effort session label from path ...\subject\date\session
sp = split(string(datpath), filesep);
if numel(sp) >= 3
    nm = sprintf('%s | %s | %s', sp(end-2), sp(end-1), sp(end));
else
    nm = char(datpath);
end
end
