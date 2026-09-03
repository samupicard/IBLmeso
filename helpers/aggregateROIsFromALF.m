function T = aggregateROIsFromALF(rootPaths, varargin)
% aggregateROIsFromALF  Aggregate per-FOV ROI TSVs into one table
%
% T = aggregateROIsFromALF(rootPaths, Name,Value,...)
%
% Required:
%   rootPaths : char or cell array of subject/date/session root paths
%
% Name-Value:
%   'columnName'    : name of column in mpciROIs.taskTunedP.tsv to use (required)
%   'alfSubfolder'  : default 'alf'
%   'alpha'         : two-tailed alpha (default 0.05)
%   'useType'       : false/true (default false)
%   'typeValue'     : integer value to keep from mpciROIs.mpciROITypes.npy (default 1)
%
% Behavior / assumptions:
% - Only loads mpciROIs.taskTunedP.tsv in each alf/FOV folder.
% - Loads mpciROIs.mlapdv_estimate.npy for positions (uses first two cols if >=2).
% - If useType true, loads mpciROIs..mpciROITypes.npy and keeps ROIs equal to typeValue.
% - h is two-tailed: h = p < (alpha/2) OR p > (1 - alpha/2).
% - Output table columns: subject, date, session, FOV, pos (1x2 double or []), p (double), h (logical).

% Parse inputs
p = inputParser;
addRequired(p,'rootPaths');
addParameter(p,'columnName',[],@ischar);
addParameter(p,'alfSubfolder','alf',@ischar);
addParameter(p,'alpha',0.05,@(x)isnumeric(x) && isscalar(x) && x>0 && x<1);
addParameter(p,'useType',false,@islogical);
addParameter(p,'typeValue',1,@(x)isnumeric(x) && isscalar(x));
parse(p,rootPaths,varargin{:});
opts = p.Results;

if isempty(opts.columnName)
    error('You must provide ''columnName'' corresponding to a header in mpciROIs.taskTunedP.tsv');
end

if ischar(rootPaths); rootPaths = {rootPaths}; end

FILENAME_P = 'mpciROIs.taskTunedP.tsv';
FILENAME_POS = 'mpciROIs.mlapdv_estimate.npy';
FILENAME_CLASS = 'mpciROIs.cellClassifier.npy';
FILENAME_TYPE = 'mpciROIs.mpciROITypes.npy';

rows = {}; % each row: subject,date,session,FOV,pos,p,h

for iRoot = 1:numel(rootPaths)
    base = rootPaths{iRoot};
    [subj,date,session] = splitPathThree(base);
    alfFolder = fullfile(base, opts.alfSubfolder);
    if ~isfolder(alfFolder), continue; end
    d = dir(alfFolder);
    isdirflag = [d.isdir];
    names = {d(isdirflag).name};
    names = names(~ismember(names, {'.','..'}));
    for k = 1:numel(names)
        FOVname = names{k};
        FOVpath = fullfile(alfFolder, FOVname);
        pfile = fullfile(FOVpath, FILENAME_P);
        if ~isfile(pfile), continue; end

        % Read TSV and select column
        try
            Tfile = readtable(pfile, 'FileType','text', 'Delimiter','\t', 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
        catch
            warning('Could not read %s in %s — skipping', FILENAME_P, FOVpath);
            continue;
        end
        if ~ismember(opts.columnName, Tfile.Properties.VariableNames)
            warning('Column ''%s'' not found in %s — skipping', opts.columnName, pfile);
            continue;
        end
        pcol = Tfile.(opts.columnName);
        pcol = double(pcol(:)); % ensure column vector double

        nROI = numel(pcol);

        % Load positions (npy)
        pos = [];
        posfile = fullfile(FOVpath, FILENAME_POS);
        if isfile(posfile)
            try
                posm = readNPY(posfile);
                % Expect Nx2 or Nx3 (rows = ROI)
                if size(posm,1) == nROI && size(posm,2) >= 2
                    pos = double(posm(:,1:2));
                elseif size(posm,1) == nROI && size(posm,2) == 1
                    pos = []; % single column not useful
                end
            catch
                warning('Could not read %s in %s', FILENAME_POS, FOVpath);
            end
        end

        % Optional classifier filter
        keepIdx = true(nROI,1);
        if opts.useType
            typefile = fullfile(FOVpath, FILENAME_TYPE);
            if ~isfile(typefile)
                warning('mpciROITypes requested but %s not found in %s — skipping type filter', FILENAME_TYPE, FOVpath);
            else
                try
                    typeVec = readNPY(typefile);
                    typeVec = double(typeVec(:));
                    if numel(typeVec) ~= nROI
                        warning('mpciROITypes length mismatch in %s — skipping type filter', FOVpath);
                    else
                        keepIdx = typeVec == opts.typeValue;
                    end
                catch
                    warning('Could not read %s in %s — skipping type filter', FILENAME_TYPE, FOVpath);
                end
            end
        end

        % Compute two-tailed h
        alpha2 = opts.alpha/2;
        hvec = (pcol < alpha2) | (pcol > (1 - alpha2));

        % Append rows for kept indices
        for r = find(keepIdx(:))'
            if isempty(pos)
                poscell = [];
            else
                poscell = pos(r,:);
            end
            rows(end+1,:) = { string(subj), string(date), string(session), string(FOVname), poscell, double(pcol(r)), logical(hvec(r)) }; %#ok<AGROW>
        end
    end
end

% Build table
if isempty(rows)
    T = table(string.empty, string.empty, string.empty, string.empty, {[]}', double.empty, false.empty, ...
        'VariableNames', {'subject','date','session','FOV','pos','p','h'});
    return;
end

T = table([rows{:,1}]', [rows{:,2}]', [rows{:,3}]', [rows{:,4}]', {rows{:,5}}', [rows{:,6}]', [rows{:,7}]', ...
    'VariableNames', {'subject','date','session','FOV','pos','p','h'});

end

%% Helper functions
function [subj,date,session] = splitPathThree(p)
p = normalizePathSeparators(p);
parts = strsplit(p, filesep);
if numel(parts) >= 3
    subj = parts{end-2};
    date = parts{end-1};
    session = parts{end};
else
    subj = ''; date = ''; session = '';
end
end

function s = normalizePathSeparators(p)
if ispc
    s = strrep(p, '/', '\');
else
    s = strrep(p, '\', '/');
end
end
