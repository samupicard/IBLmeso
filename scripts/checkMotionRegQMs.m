%% Batch motion-registration QC

%rootPath = 'Y:\Subjects\SP081\2026-09-15\001';   % drift
rootPath = 'Y:\Subjects\SP075\2025-11-10\001';   % no drift?


nPCsPlot = 4;
nPCsScore = 30;

% Find all Fall.mat files recursively
fallFiles = dir(fullfile(rootPath,'**','Fall.mat'));

fprintf('Found %d Fall.mat files\n',numel(fallFiles));


%% Preallocate outputs

driftScores = nan(numel(fallFiles),nPCsScore);
excursionScores = nan(numel(fallFiles),nPCsScore);
rTimes = nan(numel(fallFiles),nPCsScore);
lowFPowers = nan(numel(fallFiles),nPCsScore);

fallPaths = cell(numel(fallFiles),1);


%% Loop through files

for iFile = 1:numel(fallFiles)

    fallPath = fullfile(fallFiles(iFile).folder,fallFiles(iFile).name);

    fallPaths{iFile} = fallPath;

    fprintf('\n[%d/%d] %s\n', ...
        iFile,numel(fallFiles),fallPath);


    try

        %% Load ops only

        S = load(fallPath,'ops');

        if ~isfield(S,'ops')
            warning('No ops variable found in %s',fallPath);
            continue
        end

        ops = S.ops;


        %% Define GIF output path

        gifPath = fullfile(fallFiles(iFile).folder, ...
            'motionRegistrationMetrics.gif');


        %% Run QC

        [thisDrift,thisExcursion,thisrTime,thisLowFpower] = plotMotionRegistrationMetrics(ops, ...
            'nPCsPlot',nPCsPlot, ...
            'cropCentral', true, ...
            'saveGif',true, ...
            'gifPath',gifPath, ...
            'flipImages',false);


        %% Store drift scores

        nThis = min(numel(thisDrift),nPCsScore);
        driftScores(iFile,1:nThis) = thisDrift(1:nThis);
        excursionScores(iFile,1:nThis) = thisExcursion(1:nThis);
        rTimes(iFile,1:nThis) = thisrTime(1:nThis);
        lowFPowers(iFile,1:nThis) = thisLowFpower(1:nThis);


        %% Close figure to avoid accumulating windows

        %close(gcf);


    catch ME

        warning('Failed for:\n%s\n%s',fallPath,ME.message);

        continue

    end

end


%% Save aggregate results

% save(fullfile(rootPath,'motionRegistrationDriftScores.mat'), ...
%    'driftScores','fallPaths','nPCs');


%% Also make tables for easier inspection
varNames = compose('PC%d',1:nPCsScore);

driftScoreTable = array2table(driftScores, ...
    'VariableNames',varNames);
driftScoreTable.fallPath = fallPaths;
driftScoreTable = movevars(driftScoreTable,'fallPath','Before',1);

lowFPowerTable = array2table(lowFPowers, ...
    'VariableNames',varNames);
lowFPowerTable.fallPath = fallPaths;
lowFPowerTable = movevars(lowFPowerTable,'fallPath','Before',1);


%writetable(driftTable, ...
%    fullfile(rootPath,'motionRegistrationDriftScores.csv'));


fprintf('\nDone.\n');


%% Plot aggregate PC metrics across Fall.mat files
if true
    pcIdx = 1:nPCsScore;

    figure('Color','w', ...
        'Name','Aggregate motion registration metrics', ...
        'WindowStyle','normal','Position',[3300,0,500,600]);

    tl = tiledlayout(2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    % --------------------------------------------------------------
    % Drift score
    % --------------------------------------------------------------

    ax1 = nexttile;
    hold(ax1,'on');

    for iFile = 1:numel(fallPaths)
        plot(ax1,pcIdx,driftScores(iFile,:), ...
            '-o', ...
            'LineWidth',1, ...
            'MarkerSize',3);
    end

    ylabel(ax1,'Drift score');
    title(ax1,'Drift score');
    xlim(ax1,[1 nPCsScore]);
    %xticks(ax1,1:nPCsScore);
    ylim(ax1, [0,3]);
    grid(ax1,'on');


    % % --------------------------------------------------------------
    % % Excursion score
    % % --------------------------------------------------------------
    % 
    % ax2 = nexttile;
    % hold(ax2,'on');
    % 
    % for iFile = 1:numel(fallPaths)
    %     plot(ax2,pcIdx,excursionScores(iFile,:), ...
    %         '-o', ...
    %         'LineWidth',1, ...
    %         'MarkerSize',3);
    % end
    % 
    % ylabel(ax2,'Excursion score');
    % title(ax2,'Excursion score');
    % xlim(ax2,[1 nPCsScore]);
    % %xticks(ax2,1:nPCsScore);
    % grid(ax2,'on');
    % 
    % 
    % % --------------------------------------------------------------
    % % Absolute time correlation
    % % --------------------------------------------------------------
    % 
    % ax3 = nexttile;
    % hold(ax3,'on');
    % 
    % for iFile = 1:numel(fallPaths)
    %     plot(ax3,pcIdx,rTimes(iFile,:), ...
    %         '-o', ...
    %         'LineWidth',1, ...
    %         'MarkerSize',3);
    % end
    % 
    % ylabel(ax3,'|r(time)|');
    % xlabel(ax3,'PC');
    % title(ax3,'Absolute correlation with time');
    % xlim(ax3,[1 nPCsScore]);
    % %xticks(ax3,1:nPCsScore);
    % ylim(ax3,[0 1]);
    % grid(ax3,'on');

    % --------------------------------------------------------------
    % low F power fraction
    % --------------------------------------------------------------

    ax2 = nexttile;
    hold(ax2,'on');

    for iFile = 1:numel(fallPaths)
        plot(ax2,pcIdx,lowFPowers(iFile,:), ...
            '-o', ...
            'LineWidth',1, ...
            'MarkerSize',3);
    end

    ylabel(ax2,'low-f power fraction');
    title(ax2,'Low-f Power Fraction');
    xlim(ax2,[1 nPCsScore]);
    xlabel(ax2,'PC');
    ylim(ax2, [0,1]);
    %xticks(ax1,1:nPCsScore);
    grid(ax2,'on');
    

    % --------------------------------------------------------------
    % Legend labels from FOV folder
    % --------------------------------------------------------------

    legendLabels = cell(size(fallPaths));

    for iFile = 1:numel(fallPaths)

        parts = split(strrep(fallPaths{iFile},'\','/'),'/');

        fovIdx = find(startsWith(parts,'FOV_'),1,'last');

        if ~isempty(fovIdx)
            legendLabels{iFile} = parts{fovIdx};
        else
            legendLabels{iFile} = sprintf('File %d',iFile);
        end

    end

    legend(ax1,legendLabels, ...
        'Interpreter','none', ...
        'Location','eastoutside');




end

%% Plot drift scores with PCs sorted independently for each FOV

fig2 = figure( ...
    'Color','w', ...
    'Name','Sorted motion registration drift scores', ...
    'WindowStyle','normal',...
    'Position',[3233 77 602 397]);

ax = axes(fig2);
hold(ax,'on');

pcRank = 1:nPCsScore;

for iFile = 1:numel(fallPaths)

    sortedScores = sort(lowFPowers(iFile,:), ...
        'descend', ...
        'MissingPlacement','last');

    plot(ax,pcRank,sortedScores, ...
        '-', ...
        'LineWidth',1.5, ...
        'MarkerSize',3);

end

xlabel(ax,'PC rank');
ylabel(ax,'Drift score');

title(ax,'low-f power fraction per PC, sorted by magnitude');

xlim(ax,[1 nPCsScore]);
%xticks(ax,1:nPCsScore);

ylim(ax,[0,1]);
grid(ax,'on');


%legend
legendLabels = cell(size(fallPaths));

for iFile = 1:numel(fallPaths)

    parts = split(strrep(fallPaths{iFile},'\','/'),'/');

    fovIdx = find(startsWith(parts,'FOV_'),1,'last');

    if ~isempty(fovIdx)
        legendLabels{iFile} = parts{fovIdx};
    else
        legendLabels{iFile} = sprintf('File %d',iFile);
    end

end

legend(ax,legendLabels, ...
    'Interpreter','none', ...
    'Location','eastoutside');

