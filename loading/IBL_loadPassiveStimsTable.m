function passiveStimT = IBL_loadPassiveStimsTable(datpath, varargin)
% IBL_loadPassiveStimsTable Loads _ibl_passiveStims.table.csv files.
%
% passiveStimT = IBL_loadPassiveStimsTable(datpath)
%
% Searches for passive-stimulation tables in:
%
%   <datpath>\alf\task*
%
% If no task folders are found at datpath, the function searches the
% corresponding session under rootLocal.
%
% The CSV's first unnamed column is treated as row names and is not
% included as a table variable.
%
% Output columns:
%   valveOn
%   valveOff
%   toneOn
%   toneOff
%   noiseOn
%   noiseOff
%   boutIdx
%
% Example:
%
% passiveStimT = IBL_loadPassiveStimsTable( ...
%     'Y:\Subjects\SP080\2026-07-07\002');

p = inputParser;
p.addParameter('rootLocal', ...
    'C:\Users\Samuel\Documents\2PI\', ...
    @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

rootLocal = char(p.Results.rootLocal);
datpath   = char(datpath);

datpath = stripTrailingFilesep(datpath);

[subject, dateStr, session] = parseSessionPath(datpath);
datpathLocal = fullfile(rootLocal, subject, dateStr, session);

allFields = ["valveOn", "valveOff", ...
             "toneOn",  "toneOff", ...
             "noiseOn", "noiseOff"];

serverTaskFolders = dir(fullfile(datpath,      'alf', 'task*'));
localTaskFolders  = dir(fullfile(datpathLocal, 'alf', 'task*'));

taskFolders = chooseTaskSource(serverTaskFolders, localTaskFolders);

passiveStimT = array2table(zeros(0, numel(allFields)), ...
    'VariableNames', cellstr(allFields));

boutIdx = zeros(0, 1);

for iBout = 1:numel(taskFolders)
    boutpath = fullfile( ...
        taskFolders(iBout).folder, ...
        taskFolders(iBout).name);

    passiveStimFile = fullfile( ...
        boutpath, '_ibl_passiveStims.table.csv');

    if ~exist(passiveStimFile, 'file')
        continue
    end

    boutTable = readPassiveStimTable(passiveStimFile);

    boutTable = ensureColumns(boutTable, allFields);

    passiveStimT = [ ...
        passiveStimT; ...
        boutTable(:, cellstr(allFields))];

    boutIdx = [ ...
        boutIdx; ...
        iBout * ones(height(boutTable), 1)];
end

passiveStimT.boutIdx = boutIdx;

end


function datpath = stripTrailingFilesep(datpath)
while ~isempty(datpath) && ...
        (datpath(end) == '/' || datpath(end) == '\')
    datpath(end) = [];
end
end


function [subject, dateStr, session] = parseSessionPath(datpath)
pathParts = split(string(datpath), ["/", "\"]);
pathParts(pathParts == "") = [];

if numel(pathParts) < 3
    error('IBL_loadPassiveStimsTable:InvalidSessionPath', ...
        ['datpath must end with the session hierarchy ' ...
         '<subject>\<date>\<session>.']);
end

subject = char(pathParts(end-2));
dateStr = char(pathParts(end-1));
session = char(pathParts(end));
end


function taskFolders = chooseTaskSource( ...
        serverTaskFolders, localTaskFolders)

if ~isempty(serverTaskFolders)
    taskFolders = serverTaskFolders;
else
    taskFolders = localTaskFolders;
end

% Exclude non-folder matches, if any.
taskFolders = taskFolders([taskFolders.isdir]);

% Sort folders by name so that boutIdx follows task-folder order.
if ~isempty(taskFolders)
    [~, sortIdx] = sort({taskFolders.name});
    taskFolders = taskFolders(sortIdx);
end
end


function passiveStimTable = readPassiveStimTable(passiveStimFile)

% The CSV begins with an unnamed index column:
%
%   ,valveOn,valveOff,...
%   0,2621.269,2621.372,...
%
% ReadRowNames prevents that index from becoming Var1.

passiveStimTable = readtable(passiveStimFile, 'ReadRowNames',true);

% Convert any imported time columns to numeric where possible.
variableNames = passiveStimTable.Properties.VariableNames;

for iVar = 1:numel(variableNames)
    variableName = variableNames{iVar};
    values = passiveStimTable.(variableName);

    if iscell(values) || isstring(values) || ischar(values)
        numericValues = str2double(string(values));

        if all(~isnan(numericValues) | ismissing(string(values)))
            passiveStimTable.(variableName) = numericValues;
        end
    end
end
end


function inputTable = ensureColumns(inputTable, fields)
for iField = 1:numel(fields)
    field = fields(iField);

    if ~ismember(field, string(inputTable.Properties.VariableNames))
        inputTable.(field) = nan(height(inputTable), 1);
    end
end
end
