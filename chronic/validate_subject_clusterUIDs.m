function report = validate_subject_clusterUIDs(subject, varargin)
% validate_subject_clusterUIDs Validate existing mpciROIs.clusterUIDs.csv files.
%
% report = validate_subject_clusterUIDs(subject, 'Root', 'Y:\Subjects\', ...
%                                       'Dates', {}, 'ErrorOnMismatch', false)
%
% Dates: optional cellstr of date strings to filter sessions (simple substring match on paths)
% ErrorOnMismatch: if true, throws on first mismatch; otherwise collects report.

p = inputParser;
p.addParameter('Root', 'Y:\Subjects\', @(x)ischar(x)||isstring(x));
p.addParameter('Dates', {}, @(x)iscell(x)||isstring(x));
p.addParameter('ErrorOnMismatch', false, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});

root = char(p.Results.Root);
dates = cellstr(p.Results.Dates);
errorOnMismatch = p.Results.ErrorOnMismatch;

% Find all clusterUID CSVs under the subject
subjPath = fullfile(root, subject);
csvs = dir(fullfile(subjPath, '**', 'mpciROIs.clusterUIDs.csv'));

report = struct('folder', {}, 'nCSV', {}, 'nCC', {}, 'ok', {}, 'message', {});
ri = 0;

for k = 1:numel(csvs)
    folder = csvs(k).folder;
    fullCsv = fullfile(folder, csvs(k).name);

    % Optional date filter (lightweight / robust enough for most folder conventions)
    if ~isempty(dates)
        hit = false;
        for d = 1:numel(dates)
            if contains(folder, dates{d})
                hit = true; break;
            end
        end
        if ~hit
            continue
        end
    end

    ccPath = fullfile(folder, 'mpciROIs.cellClassifier.npy');
    if ~isfile(ccPath)
        msg = 'Missing mpciROIs.mpciROITypes.npy';
        if errorOnMismatch, error('ROICaT:MissingMpciROITypes', '%s (%s)', msg, folder); end
        ri = ri + 1;
        report(ri) = struct('folder', folder, 'n_clusterUIDs', NaN, 'n_mpciROIs', NaN, 'ok', false, 'message', msg);
        continue
    end

    % Read lengths
    CUIDs = read_uid_csv(fullCsv);
    %CUIDs = readmatrix(fullCsv, 'OutputType', 'string'); %THIS SKIPS EMPTY ROWS!
    nCSV = numel(CUIDs);

    cc = readNPY(ccPath);
    nCC = size(cc,1);
    if isvector(cc), nCC = numel(cc); end

    ok = (nCSV == nCC);
    msg = '';
    if ~ok
        msg = sprintf('Mismatch: clusterUIDs.csv has %d rows, mpciROITypes.npy has %d', nCSV, nCC);
        if errorOnMismatch
            error('ROICaT:CUIDLengthMismatch', '%s (folder=%s)', msg, folder);
        end
    end

    ri = ri + 1;
    report(ri) = struct('folder', folder, 'nCSV', nCSV, 'nCC', nCC, 'ok', ok, 'message', msg);
end

% Print a short summary
nChecked = numel(report);
nBad = sum(~[report.ok]);
fprintf('Checked %d folders. %d mismatches.\n', nChecked, nBad);
if nBad > 0
    fprintf('First few mismatches:\n');
    badIdx = find(~[report.ok]);
    for i = 1:min(10, numel(badIdx))
        r = report(badIdx(i));
        fprintf(' - %s | %s\n', r.folder, r.message);
    end
end
end
