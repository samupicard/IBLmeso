function T = aggregateROIsFromALF_fast(rootPaths, varargin)
%AGGREGATEROISFROMALF_FAST Aggregate ROI metadata and analysis outputs.
%
% T = aggregateROIsFromALF_fast(rootPaths, Name,Value,...)
%
% The function aggregates ROI-level data from all FOV folders below one or
% more session paths.
%
% It can operate in two modes:
%
%   Task-stat mode
%   --------------
%   Supply:
%
%       'columnName', statName
%
%   The function loads the requested columns from:
%
%       mpciROIs.taskTunedP.tsv
%       mpciROIs.taskTunedStat.tsv
%
%   Passive-movie-only mode
%   -----------------------
%   Leave 'columnName' empty:
%
%       'columnName', []
%
%   The function then uses:
%
%       mpciROIs.passiveMovieCorr.npy
%
%   to determine the number of ROIs and does not require task-tuning TSVs.
%   Task-related output columns are filled with NaN or false.
%
% Optional name-value arguments
% -----------------------------
% 'columnName'     Task statistic name. Default: [].
% 'alfSubfolder'   ALF subfolder name. Default: 'alf'.
% 'alpha'          Significance level. Default: 0.05.
% 'useChronic'     Keep only ROIs with valid chronic UIDs.
% 'useType'        Apply ROI-type filtering.
% 'typeValue'      ROI type to retain. Default: 1.
% 'onlyResponsive' Apply task-responsiveness filtering.
%
% Output columns
% --------------
% subject, date, session, FOV, cUID, ML, AP,
% stat, stat_odd, stat_even, p, h, passiveMovieCorr


p = inputParser;

addRequired(p, 'rootPaths');

addParameter( ...
    p, 'columnName', [], ...
    @(x) isempty(x) || ischar(x) || isstring(x));

addParameter(p, 'alfSubfolder', 'alf', @ischar);

addParameter( ...
    p, 'alpha', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);

addParameter(p, 'useChronic', false, @islogical);
addParameter(p, 'useType', false, @islogical);
addParameter(p, 'onlyResponsive', false, @islogical);

addParameter( ...
    p, 'typeValue', 1, ...
    @(x) isnumeric(x) && isscalar(x));

parse(p, rootPaths, varargin{:});
opts = p.Results;

opts.loadTaskStats = ~isempty(opts.columnName);

if opts.loadTaskStats
    opts.columnName = char(opts.columnName);
else
    opts.columnName = '';
end

if ischar(rootPaths) || isstring(rootPaths)
    rootPaths = cellstr(rootPaths);
end

% ----- Filenames and constants -----
FNAME_TASK_STAT    = 'mpciROIs.taskTunedStat.tsv';
FNAME_TASK_P       = 'mpciROIs.taskTunedP.tsv';
FNAME_PASSIVE_STAT = 'mpciROIs.passiveTunedStat.tsv';
FNAME_PASSIVE_P    = 'mpciROIs.passiveTunedP.tsv';
FNAME_POS   = 'mpciROIs.mlapdv_estimate.npy';
FNAME_TYPE  = 'mpciROIs.mpciROITypes.npy';
FNAME_RESP = 'mpciROIs.taskResponsiveP.tsv';
FNAME_CUID  = 'mpciROIs.clusterUIDs.csv';
FNAME_MOVIE  = 'mpciROIs.passiveMovieCorr.npy';

alpha2 = opts.alpha / 2;

% ----- Collect all FOV paths -----
allFOV = cell(0,5); % {FOVpath, subj, date, session, FOVname}
for i = 1:numel(rootPaths)
    base = rootPaths{i};
    [subj,date,session] = splitPathThree(base);
    alfDir = fullfile(base, opts.alfSubfolder);
    if ~isfolder(alfDir)
        continue;
    end

    d = dir(alfDir);
    isdirflag = [d.isdir];
    names = {d(isdirflag).name};
    names = names(~ismember(names,{'.','..'}));

    for k = 1:numel(names)
        allFOV(end+1,:) = {fullfile(alfDir,names{k}), subj, date, session, names{k}}; %#ok<AGROW>
    end
end

nFOV = size(allFOV,1);
if nFOV == 0
    T = emptyOutputTable();
    return;
end

% ----- Process FOVs in parallel if possible -----
useParallel = true;
if isempty(gcp('nocreate'))
    parpool('threads');
end

fovResults = cell(nFOV,1);

if useParallel
    parfor i = 1:nFOV
        fovResults{i} = processFOV( ...
            allFOV(i,:), opts, ...
            FNAME_TASK_STAT, FNAME_TASK_P, ...
            FNAME_PASSIVE_STAT, FNAME_PASSIVE_P, ...
            FNAME_POS, FNAME_TYPE, ...
            FNAME_CUID, FNAME_RESP, FNAME_MOVIE, alpha2);
    end
else
    for i = 1:nFOV
        fovResults{i} = processFOV( ...
            allFOV(i,:), opts, ...
            FNAME_TASK_STAT, FNAME_TASK_P, ...
            FNAME_PASSIVE_STAT, FNAME_PASSIVE_P, ...
            FNAME_POS, FNAME_TYPE, ...
            FNAME_CUID, FNAME_RESP, FNAME_MOVIE, alpha2);
    end
end




% ----- Concatenate results once -----
subjects        = {};
dates           = {};
sessions        = {};
FOVs            = {};
cUIDs           = strings(0,1);
posMs           = [];
stats           = [];
statsOdd        = [];
statsEven       = [];
pvals           = [];
hvals           = false(0,1);
passiveMovieCorr = [];

for i = 1:nFOV
    R = fovResults{i};
    if isempty(R)
        continue;
    end

    n = numel(R.passiveMovieCorr);
    subjects(end+1:end+n)   = R.subject;    %#ok<AGROW>
    dates(end+1:end+n)      = R.date;       %#ok<AGROW>
    sessions(end+1:end+n)   = R.session;    %#ok<AGROW>
    FOVs(end+1:end+n)       = R.FOV;        %#ok<AGROW>
    cUIDs(end+1:end+n)      = R.cUID;       %#ok<AGROW>
    posMs(end+1:end+n,:)    = R.pos;        %#ok<AGROW>
    stats(end+1:end+n,1)    = R.stat;       %#ok<AGROW>
    statsOdd(end+1:end+n,1) = R.stat_odd;   %#ok<AGROW>
    statsEven(end+1:end+n,1)= R.stat_even;  %#ok<AGROW>
    pvals(end+1:end+n,1)    = R.p;          %#ok<AGROW>
    hvals(end+1:end+n,1)    = R.h;          %#ok<AGROW>
    passiveMovieCorr(end+1:end+n,1) = R.passiveMovieCorr; %#ok<AGROW>
end

if isempty(passiveMovieCorr)
    T = emptyOutputTable();
else
    T = table( ...
        string(subjects(:)), ...
        string(dates(:)), ...
        string(sessions(:)), ...
        string(FOVs(:)), ...
        string(cUIDs(:)), ...
        posMs(:,1), ...
        posMs(:,2), ...
        stats(:), ...
        statsOdd(:), ...
        statsEven(:), ...
        pvals(:), ...
        hvals(:), ...
        passiveMovieCorr(:), ...
        'VariableNames', { ...
        'subject', ...
        'date', ...
        'session', ...
        'FOV', ...
        'cUID', ...
        'ML', ...
        'AP', ...
        'stat', ...
        'stat_odd', ...
        'stat_even', ...
        'p', ...
        'h', ...
        'passiveMovieCorr'});
end

end

% ===== Helper: process a single FOV =====
function R = processFOV( ...
    fovRow, opts, ...
    FNAME_TASK_STAT, FNAME_TASK_P, ...
    FNAME_PASSIVE_STAT, FNAME_PASSIVE_P, ...
    FNAME_POS, FNAME_TYPE, ...
    FNAME_CUID, FNAME_RESP, FNAME_MOVIE, alpha2)
%PROCESSFOV Load ROI data from one FOV.
%
% Supports:
%   1. Task-stat mode when columnName exists in task-tuning TSVs.
%   2. Passive-stat mode when columnName exists in passive-tuning TSVs.
%   3. Passive-movie-only mode when columnName is empty.
%
% If a requested column exists in both task and passive files, the task
% files take precedence to preserve the previous behaviour.
%
% Task p-values are percentile-style values:
%   significant when p < alpha/2 or p > 1-alpha/2.
%
% Passive p-values are conventional two-sided permutation p-values:
%   significant when p < alpha.

R = [];

FOVpath = fovRow{1};
subj    = string(fovRow{2});
date    = string(fovRow{3});
session = string(fovRow{4});
FOVname = string(fovRow{5});


%% Load passive-movie correlation, if available

movieFile = fullfile(FOVpath, FNAME_MOVIE);
passiveMovieCorr = [];

if isfile(movieFile)
    try
        passiveMovieCorr = single(readNPY(movieFile));
        passiveMovieCorr = passiveMovieCorr(:);
    catch ME
        warning( ...
            'Could not read %s (%s).', ...
            movieFile, ME.message);

        passiveMovieCorr = [];
    end
end


%% Determine ROI count and load requested statistic

if opts.loadTaskStats

    statMainName = opts.columnName;
    statOddName  = [opts.columnName '_odd'];
    statEvenName = [opts.columnName '_even'];

    taskPFile = fullfile(FOVpath, FNAME_TASK_P);
    taskStatFile = fullfile(FOVpath, FNAME_TASK_STAT);

    passivePFile = fullfile(FOVpath, FNAME_PASSIVE_P);
    passiveStatFile = fullfile(FOVpath, FNAME_PASSIVE_STAT);

    % Determine which table pair contains the requested column.
    [statSource, pfile, statfile] = chooseStatisticSource( ...
        statMainName, ...
        taskPFile, taskStatFile, ...
        passivePFile, passiveStatFile);

    if statSource == ""
        return
    end

    % Load p-values.
    try
        pOpts = detectImportOptions( ...
            pfile, ...
            'FileType', 'text', ...
            'Delimiter', '\t', ...
            'VariableNamingRule', 'preserve');

        if ~ismember(statMainName, pOpts.VariableNames)
            return
        end

        pOpts.SelectedVariableNames = {statMainName};
        Tpcol = readtable(pfile, pOpts);

        pcol = single(Tpcol.(statMainName));
        pcol = pcol(:);

    catch ME
        warning( ...
            'Could not load column %s from %s (%s).', ...
            statMainName, pfile, ME.message);
        return
    end

    nROI = numel(pcol);

    if nROI == 0
        return
    end

    % Load statistics.
    try
        sOpts = detectImportOptions( ...
            statfile, ...
            'FileType', 'text', ...
            'Delimiter', '\t', ...
            'VariableNamingRule', 'preserve');

        if ~ismember(statMainName, sOpts.VariableNames)
            return
        end

        if statSource == "task"
            % Task tables may contain cross-validation split columns.
            hasOdd = ismember(statOddName, sOpts.VariableNames);
            hasEven = ismember(statEvenName, sOpts.VariableNames);

            selectedStatVars = {statMainName};

            if hasOdd
                selectedStatVars{end+1} = statOddName;
            end

            if hasEven
                selectedStatVars{end+1} = statEvenName;
            end

        else
            % Passive tuning tables contain only the main statistic.
            hasOdd = false;
            hasEven = false;
            selectedStatVars = {statMainName};
        end

        sOpts.SelectedVariableNames = selectedStatVars;
        Tscol = readtable(statfile, sOpts);

        stat = single(Tscol.(statMainName));
        stat = stat(:);

        if hasOdd
            stat_odd = single(Tscol.(statOddName));
            stat_odd = stat_odd(:);
        else
            stat_odd = nan(nROI, 1, 'single');
        end

        if hasEven
            stat_even = single(Tscol.(statEvenName));
            stat_even = stat_even(:);
        else
            stat_even = nan(nROI, 1, 'single');
        end

    catch ME
        warning( ...
            'Could not load column %s from %s (%s).', ...
            statMainName, statfile, ME.message);
        return
    end

    if numel(stat) ~= nROI || ...
            numel(stat_odd) ~= nROI || ...
            numel(stat_even) ~= nROI

        warning( ...
            'Statistic and p-value row-count mismatch in %s.', ...
            FOVpath);
        return
    end

    % Task and passive tables use different p-value conventions.
    if statSource == "task"
        % Existing task-tuning percentile convention.
        hvec = ...
            (pcol < alpha2) | ...
            (pcol > 1 - alpha2);
    else
        % Conventional two-sided permutation p-value.
        hvec = pcol < opts.alpha;
    end

else

    % In movie-only mode, the passive-movie file defines the ROI count.
    if isempty(passiveMovieCorr)
        return
    end

    nROI = numel(passiveMovieCorr);

    stat      = nan(nROI, 1, 'single');
    stat_odd  = nan(nROI, 1, 'single');
    stat_even = nan(nROI, 1, 'single');
    pcol      = nan(nROI, 1, 'single');
    hvec      = false(nROI, 1);

end


%% Validate passive-movie vector against ROI count

if isempty(passiveMovieCorr)

    passiveMovieCorr = nan(nROI, 1, 'single');

elseif numel(passiveMovieCorr) ~= nROI

    warning( ...
        ['Passive-movie correlation size mismatch in %s: ' ...
        'expected %d ROIs, found %d. Filling with NaNs.'], ...
        FOVpath, nROI, numel(passiveMovieCorr));

    passiveMovieCorr = nan(nROI, 1, 'single');

end


%% Load chronic UIDs

uidCol = repmat("", nROI, 1);
hasValidUIDs = false;

cuidFile = fullfile(FOVpath, FNAME_CUID);

if isfile(cuidFile)
    try
        tmpUID = string(read_uid_csv(cuidFile));
        tmpUID = tmpUID(:);

        if numel(tmpUID) == nROI
            uidCol = tmpUID;
            hasValidUIDs = true;
        else
            warning( ...
                ['Size mismatch between UID and statistic data in %s; ' ...
                'filling empty UIDs for this FOV.'], ...
                FOVpath);
        end

    catch
        warning( ...
            ['Could not read UID file in %s; ' ...
            'filling empty UIDs for this FOV.'], ...
            FOVpath);
    end
end

if opts.useChronic && ~hasValidUIDs
    return
end


%% Load ROI positions

pos = nan(nROI, 2, 'single');
posfile = fullfile(FOVpath, FNAME_POS);

if isfile(posfile)
    try
        posm = readNPY(posfile);

        if size(posm, 1) == nROI && size(posm, 2) >= 2
            pos = single(posm(:, 1:2));
        end
    catch
        % Leave positions as NaN.
    end
end


%% Optional ROI-type filter

keepIdx = true(nROI, 1);

if opts.useType

    typefile = fullfile(FOVpath, FNAME_TYPE);

    if isfile(typefile)
        try
            typeVec = readNPY(typefile);
            typeVec = double(typeVec(:));

            if numel(typeVec) == nROI
                keepIdx = typeVec == opts.typeValue;
            else
                keepIdx = false(nROI, 1);
            end

        catch
            keepIdx = false(nROI, 1);
        end

    else
        keepIdx = false(nROI, 1);
    end
end


%% Optional task-responsive filter

if opts.onlyResponsive

    respKeep = true(nROI, 1);
    respFile = fullfile(FOVpath, FNAME_RESP);

    if isfile(respFile)
        try
            respTab = readtable( ...
                respFile, ...
                'FileType', 'text', ...
                'Delimiter', '\t', ...
                'VariableNamingRule', 'preserve');

            Resp = table2array(respTab);

            if size(Resp, 1) == nROI
                respKeep = any( ...
                    Resp < (0.001 / size(Resp, 2)), ...
                    2);
            else
                warning( ...
                    ['Responsive table size mismatch in %s; ' ...
                    'using all ROIs.'], ...
                    FOVpath);
            end

        catch ME
            warning( ...
                'Could not load %s (%s). Using all ROIs.', ...
                respFile, ME.message);
        end
    end

    keepIdx = keepIdx & respKeep;
end


%% Apply chronic filtering

if opts.useChronic
    keepIdx = keepIdx & uidCol ~= "";
end


%% Select retained ROIs

idxKeep = find(keepIdx);

if isempty(idxKeep)
    return
end

R.subject = repmat({subj}, numel(idxKeep), 1);
R.date = repmat({date}, numel(idxKeep), 1);
R.session = repmat({session}, numel(idxKeep), 1);
R.FOV = repmat({FOVname}, numel(idxKeep), 1);

R.cUID = uidCol(idxKeep);
R.pos = pos(idxKeep, :);

R.stat = stat(idxKeep);
R.stat_odd = stat_odd(idxKeep);
R.stat_even = stat_even(idxKeep);

R.p = pcol(idxKeep);
R.h = hvec(idxKeep);

R.passiveMovieCorr = passiveMovieCorr(idxKeep);

end


function [source, pfile, statfile] = chooseStatisticSource( ...
    columnName, ...
    taskPFile, taskStatFile, ...
    passivePFile, passiveStatFile)
%CHOOSESTATISTICSOURCE Find the table pair containing columnName.
%
% Task files are checked first to preserve the original behaviour.

source = "";
pfile = "";
statfile = "";

if tablePairContainsColumn( ...
        taskPFile, ...
        taskStatFile, ...
        columnName)

    source = "task";
    pfile = taskPFile;
    statfile = taskStatFile;
    return
end

if tablePairContainsColumn( ...
        passivePFile, ...
        passiveStatFile, ...
        columnName)

    source = "passive";
    pfile = passivePFile;
    statfile = passiveStatFile;
end

end


function tf = tablePairContainsColumn( ...
    pfile, ...
    statfile, ...
    columnName)
%TABLEPAIRCONTAINSCOLUMN Check that both tables contain a requested column.

tf = false;

if ~isfile(pfile) || ~isfile(statfile)
    return
end

try
    pOpts = detectImportOptions( ...
        pfile, ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'VariableNamingRule', 'preserve');

    if ~ismember(columnName, pOpts.VariableNames)
        return
    end

    sOpts = detectImportOptions( ...
        statfile, ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'VariableNamingRule', 'preserve');

    tf = ismember(columnName, sOpts.VariableNames);

catch
    tf = false;
end

end

% ===== Utility helpers =====
function [subj,date,session] = splitPathThree(p)
p = normalizeSep(p);
parts = strsplit(p, filesep);
if numel(parts) >= 3
    subj = parts{end-2};
    date = parts{end-1};
    session = parts{end};
else
    subj = '';
    date = '';
    session = '';
end
end

function s = normalizeSep(p)
if ispc
    s = strrep(p, '/', '\');
else
    s = strrep(p, '\', '/');
end
end

function T = emptyOutputTable()

T = table( ...
    string.empty(0,1), ...
    string.empty(0,1), ...
    string.empty(0,1), ...
    string.empty(0,1), ...
    string.empty(0,1), ...
    single.empty(0,1), ...
    single.empty(0,1), ...
    single.empty(0,1), ...
    single.empty(0,1), ...
    single.empty(0,1), ...
    single.empty(0,1), ...
    false(0,1), ...
    single.empty(0,1), ...
    'VariableNames', { ...
    'subject', ...
    'date', ...
    'session', ...
    'FOV', ...
    'cUID', ...
    'ML', ...
    'AP', ...
    'stat', ...
    'stat_odd', ...
    'stat_even', ...
    'p', ...
    'h', ...
    'passiveMovieCorr'});

end