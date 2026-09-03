function trialsT = IBL_loadTrialsTable(datpath, varargin)
% IBL_loadTrialsTable Loads _ibl_trials.table.pqt into a MATLAB table.
%
% Supports trainingChoiceWorld and biasedChoiceWorld data.
% Can also load bpod-time trials if requested and available locally.

p = inputParser;
p.addParameter('rootLocal', 'C:\Users\Samuel\Documents\2PI\', @ischar);
p.addParameter('sync', 'timeline', @ischar);
p.parse(varargin{:});

rootLocal = p.Results.rootLocal;
syncFlag  = p.Results.sync;

if strcmp(datpath(end), filesep)
    datpath = datpath(1:end-1);
end

[pathParts, subject, dateStr, session] = parseSessionPath(datpath);
datpathLocal = fullfile(rootLocal, subject, dateStr, session);

allFields = ["intervals_0", "intervals_1", "goCue_times", ...
    "response_times", "choice", "stimOn_times", "contrastLeft", ...
    "contrastRight", "feedback_times", "feedbackType", "rewardVolume", ...
    "probabilityLeft", "firstMovement_times", "quiescenceOn_times", ...
    "choiceMovement_times"];

serverTaskFolders = dir(fullfile(datpath,      'alf', 'task*'));
localTaskFolders  = dir(fullfile(datpathLocal, 'alf', 'task*'));

[taskFolders, trialsFilename, tableFields] = chooseTaskSource( ...
    serverTaskFolders, localTaskFolders, syncFlag, allFields);

trialsT = array2table(zeros(0, numel(allFields)), ...
    'VariableNames', cellstr(allFields));

boutIdx = [];

for iBout = 1:numel(taskFolders)
    boutpath = fullfile(taskFolders(iBout).folder, taskFolders(iBout).name);
    trialFile = fullfile(boutpath, trialsFilename);

    if ~exist(trialFile, 'file')
        continue
    end

    trialsTable = readTrialsTable(trialFile, tableFields);

    trialsTable = addQuiescenceOnTimes(trialsTable, boutpath, localTaskFolders, iBout);

    trialsTable = addChoiceMovementTimes(trialsTable, boutpath);

    trialsTable = ensureColumns(trialsTable, allFields);

    boutIdx = [boutIdx; iBout * ones(height(trialsTable), 1)];
    trialsT = [trialsT; trialsTable(:, cellstr(allFields))];
end

trialsT.contrastDiff = fillmissingWithZero(trialsT.contrastRight) - ...
    fillmissingWithZero(trialsT.contrastLeft);

trialsT.boutIdx = boutIdx;

end

function [pathParts, subject, dateStr, session] = parseSessionPath(datpath)
pathParts = split(datpath, filesep);
subject  = pathParts{end-2};
dateStr  = pathParts{end-1};
session  = pathParts{end};
end

function [taskFolders, trialsFilename, tableFields] = chooseTaskSource( ...
    serverTaskFolders, localTaskFolders, syncFlag, allFields)

if isempty(serverTaskFolders)
    taskFolders = localTaskFolders;

    if strcmp(syncFlag, 'bpod')
        trialsFilename = '_ibl_trials.table_bpod.npy';
    else
        trialsFilename = '_ibl_trials.table.pqt';
    end

    tableFields = allFields(1:end-1);
else
    taskFolders = serverTaskFolders;
    trialsFilename = '_ibl_trials.table.pqt';
    tableFields = allFields;
end
end

function trialsTable = readTrialsTable(trialFile, tableFields)
try
    trialsTable = parquetread(trialFile);
catch
    trialsTable = array2table(readNPY(trialFile), ...
        'VariableNames', cellstr(tableFields));
end
end

function trialsTable = addQuiescenceOnTimes( ...
    trialsTable, boutpath, localTaskFolders, iBout)

if any(strcmp('quiescenceOn_times', trialsTable.Properties.VariableNames))
    return
end

nTrials = height(trialsTable);
quiescenceOn_times = nan(nTrials, 1);

quiescencePeriodFile = fullfile(boutpath, '_ibl_trials.quiescencePeriod.npy');

if exist(quiescencePeriodFile, 'file')
    quiescencePeriod = readNPY(quiescencePeriodFile);
    quiescenceOn_times = trialsTable.goCue_times - quiescencePeriod;

elseif ~isempty(localTaskFolders) && iBout <= numel(localTaskFolders)
    localBoutPath = fullfile(localTaskFolders(iBout).folder, ...
        localTaskFolders(iBout).name);

    quiescenceFile = fullfile(localBoutPath, ...
        '_ibl_trials.quiescenceOn_times_bpod.npy');
    goCueFile = fullfile(localBoutPath, ...
        '_ibl_trials.goCue_times_bpod.npy');

    if exist(quiescenceFile, 'file') && exist(goCueFile, 'file')
        quiescenceOn_bpod = readNPY(quiescenceFile);
        goCue_bpod = readNPY(goCueFile);

        quiescencePeriod = goCue_bpod - quiescenceOn_bpod;
        quiescenceOn_times = trialsTable.goCue_times - quiescencePeriod;
    end
end

trialsTable = addvars(trialsTable, quiescenceOn_times);
end

function trialsTable = addChoiceMovementTimes(trialsTable, boutpath)
if any(strcmp('choiceMovement_times', trialsTable.Properties.VariableNames))
    return
end

wheelFile = fullfile(boutpath, '_ibl_wheelMoves.intervals.npy');

if ~exist(wheelFile, 'file')
    return
end

wheelMovIntervals = readNPY(wheelFile);
choiceMovement_times = nan(height(trialsTable), 1);

for iTrial = 1:height(trialsTable)
    iNearestWheelMov = find( ...
        wheelMovIntervals(:, 2) > trialsTable.response_times(iTrial), 1);

    if ~isempty(iNearestWheelMov)
        choiceMovement_times(iTrial) = wheelMovIntervals(iNearestWheelMov, 1);
    end
end

trialsTable = addvars(trialsTable, choiceMovement_times);
end

function trialsTable = ensureColumns(trialsTable, fields)
for iField = 1:numel(fields)
    field = fields(iField);

    if ~any(strcmp(field, trialsTable.Properties.VariableNames))
        trialsTable.(field) = nan(height(trialsTable), 1);
    end
end
end

function x = fillmissingWithZero(x)
x(isnan(x)) = 0;
end