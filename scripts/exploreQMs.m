%% prep
%datpath = 'Y:\Subjects\SP081\2026-09-10\001';
datpath = 'Y:\Subjects\SP076\2025-11-10\001';

%load
splitPaths = split(datpath,filesep);
subj = splitPaths{end-2};
date = splitPaths{end-1};
sess = splitPaths{end};

fprintf('%s: loading ',datpath)
Fall = IBL_loadMesoData(subj,date,sess,'fast',true);
F = Fall.F';
Fneu = Fall.Fneu';
times = Fall.time;

%% do the QC

fprintf('\nComputing QC metrics..');
neuralQMs = get_neuralQMs(F,Fneu,times);
passTable = threshold_neuralQMs(neuralQMs);
fprintf('. Done!\n');

%% Find ROIs failing at least one QC

failMask = ~all(passTable{:,:}, 2) & Fall.iscell;
failROIs = find(failMask);

maxRowsPerFig = 10;
nFail = numel(failROIs);
nFigs = ceil(nFail / maxRowsPerFig);
nFigs_toplot = min(nFigs,10);

fprintf('%d / %d ROIs fail at least one QC.\n', ...
    nFail, height(passTable));

%% plot
for iFig = 1:nFigs_toplot

    idx = (iFig-1)*maxRowsPerFig + 1 : ...
        min(iFig*maxRowsPerFig, nFail);

    theseROIs = failROIs(idx);

    figure('Name', sprintf('Failed neural QC %d/%d', iFig, nFigs));

    tl = tiledlayout(numel(theseROIs), 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax = gobjects(numel(theseROIs),1);

    for i = 1:numel(theseROIs)

        iROI = theseROIs(i);

        ax(i) = nexttile;

        plot(times, F(:,iROI), ...
            'DisplayName', 'F');
        hold on

        plot(times, Fneu(:,iROI), ...
            'DisplayName', 'Fneu');

        % Find names of failed QMs
        failedQMs = passTable.Properties.VariableNames( ...
            ~passTable{iROI,:});

        failedStr = strjoin(failedQMs, ', ');

        title(sprintf('ROI %d | fails: %s', ...
            iROI, failedStr), ...
            'Interpreter', 'none');

        %ylabel('Fluorescence');

        if i == numel(theseROIs)
            xlabel('Time (s)');
        else
            ax(i).XTickLabel = [];
        end

        if i == 1
            legend('Location', 'best');
        end

        box off

    end

    title(tl, sprintf( ...
        'ROIs failing neural QC (%d-%d of %d)', ...
        idx(1), idx(end), nFail));

    linkaxes(ax,'x');
    xlim(ax(1),[100,700]);

end

%% make diagnostic table

% Convert neuralQMs to table
qmTable = struct2table(neuralQMs);

% Keep failing ROIs only
failedQMTable = qmTable(failROIs,:);

% Add original ROI index
failedQMTable = addvars(failedQMTable, failROIs, ...
    'Before', 1, 'NewVariableNames', 'roiIdx');

% Add a readable list of failed QMs
failedNames = strings(nFail,1);

for i = 1:nFail
    failedNames(i) = strjoin( ...
        passTable.Properties.VariableNames(~passTable{failROIs(i),:}), ...
        ', ');
end

failedQMTable = addvars(failedQMTable, failedNames, ...
    'After', 'roiIdx', 'NewVariableNames', 'failedQMs');

disp(failedQMTable(1:min(100,length(failROIs)),:));