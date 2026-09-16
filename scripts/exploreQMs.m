%% prep

%datpath = 'Y:\Subjects\SP081\2026-09-10\001'; %multi-depth, some saturated ROIs
datpath = 'Y:\Subjects\SP076\2025-11-10\001'; %some low SNR FOVs


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

    % Select top and bottom pair, the pair around the threshold, and the rest evenly spaced in rank
    nThis = min(nExamples,numel(rankedROIs));
    nRanked = numel(rankedROIs);

    if nThis <= 4

        exampleRanks = unique(round(linspace(1,nRanked,nThis)));

    else

        rankedPass = passTable{rankedROIs,metricName};

        % ROIs that pass every QC metric
        allPass = all(passTable{:,:},2);
        rankedAllPass = allPass(rankedROIs);
        allPassRanks = find(rankedAllPass);

        % ROIs that fail ONLY this metric
        otherMetrics = setdiff(metricNames,metricName,'stable');

        failOnlyThisMetric = ...
            ~passTable{:,metricName} & ...
            all(passTable{:,otherMetrics},2);

        rankedFailOnly = failOnlyThisMetric(rankedROIs);
        failOnlyRanks = find(rankedFailOnly);

        % Find pass/fail boundary in the good -> bad ranking
        passRanks = find(rankedPass);
        failRanks = find(~rankedPass);

        canSeparate = ...
            ~isempty(passRanks) && ...
            ~isempty(failRanks) && ...
            passRanks(end) < failRanks(1);

        if canSeparate
            if numel(allPassRanks) >= 2
                bestRanks = allPassRanks(1:2);
            else
                bestRanks = [1 2];
            end

            if ~isempty(allPassRanks)
                closestPassRank = allPassRanks(end);
            else
                closestPassRank = passRanks(end);
            end

            if ~isempty(failOnlyRanks)
                closestFailRank = failOnlyRanks(1);
                worstFailOnlyRank = failOnlyRanks(end);
            else
                closestFailRank = failRanks(1);
                worstFailOnlyRank = [];
            end

            worstRanks = [nRanked-1 nRanked];

            reservedRanks = [ ...
                bestRanks(:)', ...
                closestPassRank, ...
                closestFailRank, ...
                worstFailOnlyRank(:)', ...
                worstRanks(:)'];

            reservedRanks = unique(reservedRanks);

            % -------------------------------------------------------------
            % Fill remaining slots
            % Prefer:
            %   - ROIs passing ALL metrics
            %   - ROIs failing ONLY this metric
            % -------------------------------------------------------------

            nRemaining = nThis - numel(reservedRanks);

            if nRemaining > 0

                preferredRanks = find(rankedAllPass | rankedFailOnly);
                preferredRanks = setdiff(preferredRanks,reservedRanks);

                if numel(preferredRanks) >= nRemaining

                    idx = round(linspace(1,numel(preferredRanks),nRemaining));
                    middleRanks = preferredRanks(idx);

                else

                    % Use all preferred examples first
                    middleRanks = preferredRanks;

                    nStillNeeded = nRemaining - numel(middleRanks);

                    % Then fill any remaining slots from the full rank range
                    candidateRanks = setdiff( ...
                        1:nRanked, ...
                        [reservedRanks middleRanks]);

                    if nStillNeeded > 0 && ~isempty(candidateRanks)

                        idx = round(linspace(1,numel(candidateRanks), ...
                            min(nStillNeeded,numel(candidateRanks))));

                        middleRanks = [middleRanks(:)' candidateRanks(idx(:))'];

                    end

                end

            else
                middleRanks = [];
            end

            % Final display order always remains good -> bad
            exampleRanks = sort(unique([reservedRanks(:)' middleRanks(:)']));

        else

            % Fallback if there is no clean pass/fail separation
            edgeRanks = [1 2 nRanked-1 nRanked];

            nRemaining = nThis - numel(unique(edgeRanks));
            candidateRanks = 3:nRanked-2;

            if nRemaining > 0 && ~isempty(candidateRanks)

                idx = round(linspace(1,numel(candidateRanks), ...
                    min(nRemaining,numel(candidateRanks))));

                middleRanks = candidateRanks(idx);

            else
                middleRanks = [];
            end

            exampleRanks = sort(unique([edgeRanks middleRanks]));

        end

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
        else
            titleColor = [0.8 0 0];
        end

        otherMetrics = setdiff(metricNames,metricName,'stable');

        failsOnlyThis = ...
            ~passTable{iROI,metricName} && ...
            all(passTable{iROI,otherMetrics});

        if failsOnlyThis
            failStr = 'FAIL ONLY';
        elseif passTable{iROI,metricName}
            failStr = 'PASS';
        else
            failStr = 'FAIL';
        end

        if isCell(iROI)
            cellStr = 'CELL';
        else
            cellStr = 'NON-CELL';
        end

        title(sprintf( ...
            'ROI %d | %s | %s = %.4g | %s', ...
            iROI,cellStr,metricName,metricVal,failStr), ...
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

    nFailMetric = sum(~passTable{:,metricName});
    nTotalMetric = height(passTable);

    title(tl,sprintf( ...
        '%s: examples ranked from good to bad | %d/%d fail (%.1f%%)', ...
        metricName,nFailMetric,nTotalMetric, ...
        100*nFailMetric/nTotalMetric), ...
        'Interpreter','none');

    linkaxes(ax,'x');
    xlim(ax(1),[3000,3600]);

end

%% Plot distributions of all QC metrics

metricNames = passTable.Properties.VariableNames;
nMetrics = numel(metricNames);

nCols = 4;
nRows = ceil(nMetrics/nCols);

figure('Name','Neural QC metric distributions');

tl = tiledlayout(nRows,nCols, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for iMetric = 1:nMetrics

    metricName = metricNames{iMetric};
    vals = qmTable.(metricName);

    ax = nexttile;
    hold(ax,'on');

    % Separate cells and non-cells
    cellVals = vals(isCell & isfinite(vals));
    nonCellVals = vals(~isCell & isfinite(vals));

    % Use common bin edges for both groups
    allVals = vals(isfinite(vals));

    if isempty(allVals)
        title(metricName,'Interpreter','none');
        continue
    end

    % Robust plotting range to avoid a few extreme outliers dominating
    xLim = prctile(allVals,[0.05 99.9]);

    % Make sure the QC threshold is included in the plotting range
    threshold = thresholds.(metricName).value;
    xLim(1) = min(xLim(1),threshold);
    xLim(2) = max(xLim(2),threshold);

    if xLim(1) == xLim(2)
        xLim = [min(allVals) max(allVals)];
    end

    if xLim(1) == xLim(2)
        xLim = xLim + [-0.5 0.5];
    end

    edges = linspace(xLim(1),xLim(2),51);

    histogram(ax,cellVals,edges, ...
        'Normalization','probability', ...
        'DisplayStyle','stairs', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Cell (n=%d)',length(cellVals)));

    histogram(ax,nonCellVals,edges, ...
        'Normalization','probability', ...
        'DisplayStyle','stairs', ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('Non-cell (n=%d)',length(nonCellVals)));

    % QC threshold
    threshold = thresholds.(metricName).value;
    passIf = thresholds.(metricName).passIf;

    if ismember(passIf, {'<','<='})
        labelSide = 'left';
    else
        labelSide = 'right';
    end

    xline(ax,threshold,'--', ...
        sprintf('%s %.3g', ...
        thresholds.(metricName).passIf,threshold), ...
        'LineWidth',1.2, ...
        'LabelVerticalAlignment','top',...
        'LabelHorizontalAlignment',labelSide, ...
        'LabelOrientation','horizontal',...
        'DisplayName','threshold');

    xlim(ax,xLim);

    title(ax,metricName,'Interpreter','none');
    ylabel(ax,'Fraction of ROIs');

    box(ax,'off');

    if iMetric == 1
        legend(ax,'Location','best');
    end

end

title(tl,'Neural QC metric distributions');

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
        xlim(ax(1),[3000,3600]);

    end
end