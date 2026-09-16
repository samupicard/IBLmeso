%% prep

%datpath = 'Y:\Subjects\SP081\2026-09-10\001';
datpath = 'Y:\Subjects\SP076\2025-11-10\001';

% Load
splitPaths = split(datpath,filesep);
subj = splitPaths{end-2};
date = splitPaths{end-1};
sess = splitPaths{end};

fprintf('%s: loading ',datpath)
Fall = IBL_loadMesoData(subj,date,sess,'fast',true);

F = Fall.F';
Fneu = Fall.Fneu';
times = Fall.time;

%% get the QC

fprintf('\nComputing QC metrics..');
neuralQMs = get_neuralQMs(F,Fneu,times);
fprintf('. Done!\n');

%% threshold the metrics

fprintf('Thresholding QC metrics..');
passTable = threshold_neuralQMs(neuralQMs);
fprintf('. Done!\n');

%% classify ROIs

nROIs = size(F,2);

assert(height(passTable) == nROIs, ...
    'passTable has %d rows but fluorescence data contain %d ROIs.', ...
    height(passTable),nROIs);

% Force cell classification to be a column vector
isCell = logical(Fall.iscell(:));

assert(numel(isCell) == nROIs, ...
    'Fall.iscell has %d elements but fluorescence data contain %d ROIs.', ...
    numel(isCell),nROIs);

% True if an ROI passes every QC metric
passAllQC = all(passTable{:,:},2);

% Separate cells and non-cells
cellMask = isCell;
nonCellMask = ~isCell;

% Separate pass/fail within each class
passCellMask = cellMask & passAllQC;
failCellMask = cellMask & ~passAllQC;

passNonCellMask = nonCellMask & passAllQC;
failNonCellMask = nonCellMask & ~passAllQC;

% We only care about cells failing at least one QC
failCells = find(failCellMask);

nCells = sum(cellMask);
nFail = numel(failCells);

fprintf('%d / %d cells fail at least one QC (%.1f%%).\n', ...
    nFail,nCells,100*nFail/nCells);


%% make diagnostic table for failed cells

qmTable = struct2table(neuralQMs);

% Keep only cells failing at least one QC
failedQMTable = qmTable(failCells,:);

% Add original ROI index
failedQMTable = addvars(failedQMTable,failCells, ...
    'Before',1,'NewVariableNames','roiIdx');

% Add readable list of failed QMs
failedNames = strings(nFail,1);

for i = 1:nFail

    iROI = failCells(i);

    failedNames(i) = strjoin( ...
        passTable.Properties.VariableNames(~passTable{iROI,:}), ...
        ', ');

end

failedQMTable = addvars(failedQMTable,failedNames, ...
    'After','roiIdx','NewVariableNames','failedQMs');

disp(failedQMTable(1:min(100,nFail),:));

%% Plot example traces spanning good-to-bad values for each QC metric

qmTable = struct2table(neuralQMs);
thresholds = get_neuralQMThresholds();

metricNames = passTable.Properties.VariableNames;
isCell = logical(Fall.iscell(:));

nExamples = 10;

for iMetric = 1:numel(metricNames)

    metricName = metricNames{iMetric};
    metricVals = qmTable.(metricName);

    % Determine good -> bad direction from threshold definition
    passIf = thresholds.(metricName).passIf;

    switch passIf
        case {'<','<='}
            sortDirection = 'ascend';    % low = good

        case {'>','>='}
            sortDirection = 'descend';   % high = good

        otherwise
            error('Unknown passIf "%s".',passIf);
    end

    % Include all ROIs with finite metric values.
    validMask = isfinite(metricVals);

    % Avoid trivial secondary failures caused by very low variance
    if ~strcmp(metricName,'var')
        validMask = validMask & passTable.var;
    end

    validROIs = find(validMask);
    vals = metricVals(validROIs);

    if isempty(validROIs)
        warning('No valid values for metric %s.',metricName);
        continue
    end

    % Rank from good -> bad
    [~,order] = sort(vals,sortDirection);
    rankedROIs = validROIs(order);

    % Select top and bottom pair, and the rest evenly spaced in rank
    nThis = min(nExamples,numel(rankedROIs));
    if nThis <= 4
        exampleRanks = unique(round(linspace(1,numel(rankedROIs),nThis)));
    else
        nMiddle = nThis - 4;
        edgeRanks = [1 2 numel(rankedROIs)-1 numel(rankedROIs)];
        middleRanks = round(linspace(3, ...
            numel(rankedROIs)-2, nMiddle));
        exampleRanks = unique([1 2 middleRanks ...
            numel(rankedROIs)-1 numel(rankedROIs)]);
    end
    exampleROIs = rankedROIs(exampleRanks);

    figure('Name',sprintf('QC examples: %s',metricName));

    tl = tiledlayout(numel(exampleROIs),1, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    ax = gobjects(numel(exampleROIs),1);

    for i = 1:numel(exampleROIs)

        iROI = exampleROIs(i);

        ax(i) = nexttile;

        plot(times,F(:,iROI),'DisplayName','F');
        hold on
        plot(times,Fneu(:,iROI),'DisplayName','Fneu');

        metricVal = metricVals(iROI);
        passed = passTable{iROI,metricName};

        if passed
            titleColor = [0 0.5 0];
            passStr = 'PASS';
        else
            titleColor = [0.8 0 0];
            passStr = 'FAIL';
        end

        if isCell(iROI)
            cellStr = 'CELL';
        else
            cellStr = 'NON-CELL';
        end

        title(sprintf( ...
            'ROI %d | %s | %s = %.4g | %s', ...
            iROI,cellStr,metricName,metricVal,passStr), ...
            'Interpreter','none', ...
            'Color',titleColor);

        if i == numel(exampleROIs)
            xlabel('Time (s)');
        else
            ax(i).XTickLabel = [];
        end

        if i == 1
            legend('Location','best');
        end

        box off
    end

    title(tl,sprintf('%s: examples ranked from good to bad',metricName), ...
        'Interpreter','none');

    linkaxes(ax,'x');
    xlim(ax(1),[100,700]);

end

%% plot failed cells
if false
    maxRowsPerFig = 10;

    nFigs = ceil(nFail/maxRowsPerFig);
    nFigs_toplot = min(nFigs,10);

    for iFig = 1:nFigs_toplot

        idx = (iFig-1)*maxRowsPerFig + 1 : ...
            min(iFig*maxRowsPerFig,nFail);

        theseROIs = failCells(idx);

        figure('Name',sprintf('Failed neural QC %d/%d',iFig,nFigs));

        tl = tiledlayout(numel(theseROIs),1, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        ax = gobjects(numel(theseROIs),1);

        for i = 1:numel(theseROIs)

            iROI = theseROIs(i);

            ax(i) = nexttile;

            plot(times,F(:,iROI), ...
                'DisplayName','F');
            hold on

            plot(times,Fneu(:,iROI), ...
                'DisplayName','Fneu');

            % Find names of failed QMs
            failedQMs = passTable.Properties.VariableNames( ...
                ~passTable{iROI,:});

            failedStr = strjoin(failedQMs,', ');

            title(sprintf('ROI %d | fails: %s', ...
                iROI,failedStr), ...
                'Interpreter','none');

            if i == numel(theseROIs)
                xlabel('Time (s)');
            else
                ax(i).XTickLabel = [];
            end

            if i == 1
                legend('Location','best');
            end

            box off

        end

        title(tl,sprintf( ...
            'Cells failing neural QC (%d-%d of %d)', ...
            idx(1),idx(end),nFail));

        linkaxes(ax,'x');
        xlim(ax(1),[100,700]);

    end
end