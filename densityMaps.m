
%subjLearners = {'SP035','SP037','SP044','SP046','SP054','SP058','SP060','SP061','SP063','SP066','SP067','SP072','SP075','SP076'};
%subjLearners = {'SP035','SP037','SP044','SP046','SP054','SP058','SP060','SP061'};

saveflag = true;
loadParquetIfExists = true;   % use saved allROIs_*.parquet when available
forceReaggregate = true;     % set true to ignore saved parquet files

load('canonicalSessions.mat');

%sPaths = IBL_listSessionPaths('root','Y:\Subjects','protocol',{'trainingChoiceWorld'},'mpci',true);
% sPaths = IBL_listSessionPaths('root','Y:\Subjects',...
%     'protocol',{'trainingChoiceWorld','trainingPhaseChoiceWorld','biasedChoiceWorld'},'mpci',true,...
%     'fullContrastSet',true,'trials',400,'perfOnEasy',0.8);

%sPaths = filter_session_paths_by_subject(sPaths,subjLearners);
%sPaths = filter_session_paths_by_subject(sPaths,{'SP075'}); sPaths = sPaths(1:8);
%sPaths = filter_session_paths_by_subject(sPaths,{'SP072'}); sPaths = sPaths(1:8);

%svPath = 'Y:\Subjects\results\taskVariableTuning_SingleCells\ALL';
%svPath = 'C:\Users\Samuel\Desktop\DataClub_2026-04\results\bCW';
%svPath = 'C:\Users\Samuel\Documents\2PI\mesoscope_active\analysis\PETH\densityMaps\bCW_SP072_2sessions';
%svPath = 'C:\Users\Samuel\Documents\2PI\mesoscope_active\analysis\PETH\densityMaps\passiveMovie';
svPath = 'C:\Users\Samuel\Documents\MATLAB\Code\IBLmeso\results';

% statNames = {...
%     'ccu_stimOn_0to400_stimSide100',...
%     'ccu_choiceMovement_-400to0_choice',...
%     'ccu_feedback_0to400_feedbackType'};

%% generate statNames by base/time window
%twin_ev = {[-0.6,-0.4],[-0.4,-0.2],[-0.2,0],[0,0.2],[0.2,0.4],[0.4,0.6]};
twin_ev = {[-0.3,-0.1],[-0.1,0.1],[0.1,0.3],[0.3,0.5],[0.5,0.7]};

bases = {
    {'ccMeanDiff','goCue','stimSide'}
    {'ccMeanDiff','goCue','stimSide100'}
    {'ccMeanDiff','choiceMovement','choice'}
    {'ccMeanDiff','feedback','feedbackType'}
    ...%{'ccMean','firstMovement','movement'}
    };

baseLabels = { ...
    'STIM SIDE', ...
    'STIM SIDE 100', ...
    'CHOICE', ...
    'FEEDBACK'};%, ...
%'MOVEMENT'};

%bases = {{'ccMeanDiff','goCue','stimSide100'}};
%baseLabels = {'STIM SIDE 100'};

%cmaps_frac_perBase = {[0.05,0.3],[0.05,0.3],[0.1,0.6]};
%cmaps_meds_perBase = {[0.45,0.55],[0.48,0.52],[0.45,0.55]};

%cmaps_frac_perBase = {[0.05,0.2],[0.05,0.2],[0.05,0.4]};
%cmaps_meds_perBase = {[-2,2],[-2,2],[-3,3]};

cmaps_frac_perBase = {[0.05,0.2],[0.05,0.2],[0.05,0.2],[0.05,0.4]};
cmaps_meds_perBase = {[-2,2],[-2,2],[-2,2],[-3,3]};

%cmaps_frac_perBase = {[0.05,0.2],[0.05,0.2],[0.05,0.4],[0.05,0.4]};
%cmaps_meds_perBase = {[-2,2],[-2,2],[-3,3],[-3,3]};

%% make one figure per stat/base
if true
    for b = 1:numel(bases)

        bparts = bases{b};
        saveBase = sprintf('%s_%s_%s', bparts{1}, bparts{2}, bparts{3});
        outBase = ['ROIdensity_timecourse_',saveBase];

        fig = figure( ...
            'Position',[1500,-200,900,700], ...
            'Name',outBase, ...
            'Color','k');

        nTime = numel(twin_ev);
        nCols = nTime + 1;

        tl = tiledlayout(3, nCols,...
            'TileSpacing','compact', ...
            'Padding','compact');
        %set(gcf,'Renderer','painters');

        for i = 1:numel(twin_ev)

            t = twin_ev{i} * 1000;
            t_mid = round(mean(t));

            statName = sprintf('%s_%s_%dto%d_%s', bparts{1}, bparts{2}, t(1), t(2), bparts{3});

            % make/load big table with all ROIs
            parquetPath = fullfile(svPath, ['allROIs_', statName, '.parquet']);

            if loadParquetIfExists && ~forceReaggregate && isfile(parquetPath)

                fprintf('Loading saved ROI table: %s\n', parquetPath);
                T = parquetread(parquetPath);

            else

                fprintf('Aggregating ROI table from ALF: %s\n', statName);
                T = aggregateROIsFromALF_fast( ...
                    sPaths, ...
                    'columnName', statName, ...
                    'useType', true,...
                    'onlyResponsive',false);

                if saveflag
                    if ~exist(svPath, 'dir')
                        mkdir(svPath);
                    end

                    fprintf('Saving ROI table: %s\n', parquetPath);
                    parquetwrite(parquetPath, T);
                end

            end

            %reduce table to one entry per cUID
            T_red = reduceROITableByCUID(T);

            %subsample to approach spatially uniform
            T_red_unif = spatialSubsampleROIs(T_red);

            % compute density maps
            D = IBL_computeROIdensityMaps(T);

            D.v_clim_sess = [0 10];

            % column label
            %colLabel = sprintf('%d to %d ms', t(1), t(2));
            colLabel = sprintf('%d ms',t_mid);

            col = i;

            % plot overall
            nexttile(col)
            IBL_plotDensityPanel(D, D.ratio, cmaps_frac_perBase{b});
            if i == 1
                ylabel(gca,'tot frac. sig.','Color','w','FontWeight','bold');
            end
            title(colLabel, 'Color','w');

            % plot high
            nexttile(col + nCols)
            IBL_plotDensityPanel(D, D.ratio_hi, cmaps_frac_perBase{b});
            if i == 1
                ylabel(gca,'hi frac. sig.','Color','w','FontWeight','bold');
            end

            % plot low
            nexttile(col + 2*nCols)
            IBL_plotDensityPanel(D, D.ratio_lo, cmaps_frac_perBase{b});
            if i == 1
                ylabel(gca,'lo frac. sig.','Color','w','FontWeight','bold');
            end

            % plot median stat
            % qidx = find(round(D.qvec,2) == 0.5);
            % statMed = squeeze(D.stat_quantiles(:,:,qidx));
            % nexttile(col + 3*nCols)
            % IBL_plotStatPanel(D, statMed, cmaps_meds_perBase{b})
            % if i == 1
            %     ylabel(gca,'median \Delta','Color','w','FontWeight','bold','Interpreter','tex');
            % end

            % plot scatterplot / summary map of individual neuron stats
            figMap = figure( ...
                'Position',[1500,-200,700,700], ...
                'Name',['ROISummaryMap_', statName], ...
                'Color','k');

            axMap = axes(figMap, 'Color','k');

            IBL_plotROISummaryMap(T_red_unif, 'ax', axMap, ...
                'alpha2',0.025,...
                'plotNonsigAsDots', true);

            title(axMap, sprintf('%s | %s', baseLabels{b}, colLabel), ...
                'Color','w', ...
                'Interpreter','none');

            if true
                if ~exist(svPath, 'dir')
                    mkdir(svPath);
                end

                mapOut = fullfile(svPath, ...
                    sprintf('ROISummaryMap_%s_%dms.png', saveBase, t_mid));

                exportgraphics(figMap, mapOut, ...
                    'Resolution', 300, ...
                    'BackgroundColor', 'current');
            end

            close(figMap);

        end

        % right-hand-side HSV legend, spanning visual space of the last column
        axLeg = nexttile(2*nCols);
        plotHSVLegend( ...
            axLeg, ...
            cmaps_frac_perBase{b}, ...
            D.v_clim_sess, ...
            'frac. sig.');

        % axLeg2 = nexttile(4*nCols);
        % plotStatLegend( ...
        %     axLeg2, ...
        %     cmaps_meds_perBase{b}, ...              % clim (stat range)
        %     D.v_clim_sess, ...      % brightness range
        %     'median \Delta');

        sgtitle(baseLabels{b},'Color','w','Interpreter','none');

        if saveflag
            exportgraphics(fig, fullfile(svPath,[outBase '.png']), ...
                'Resolution', 300, ...
                'BackgroundColor', 'current');

            % savefig(fig, fullfile(svPath,[outBase '.fig']));
            % set(fig, 'Renderer', 'painters');
            % exportgraphics(fig, fullfile(svPath,[outBase '.svg']), ...
            %     'ContentType', 'vector', ...
            %     'BackgroundColor', 'current');
            %print(fig, [outBase '.svg'], '-dsvg', '-painters');
        end
    end

    % %% generate statNames
    % twin_ev = {[-0.45,-0.3],[-0.3,-0.15],[-0.15,0],[0,0.15],[0.15,0.3],[0.3,0.45]};
    % bases = { ...
    %     'ccu_stimOn_stimSide', ...
    %     'ccu_choiceMovement_choice', ...
    %     'ccu_feedback_feedbackType'};
    % statNames = {};
    % for b = 1:numel(bases)
    %     for i = 1:numel(twin_ev)
    %         t = twin_ev{i} * 1000; % convert to ms
    %         statNames{end+1} = sprintf('%s_%dto%d', bases{b}, t(1), t(2));
    %     end
    % end
    %
    % %% make and save density maps
    % for i = 1:length(statNames)
    %
    %     statName = statNames{i};
    %
    %     %make big table with all ROIs
    %     T = aggregateROIsFromALF_fast(sPaths,'columnName',statName,'useType',true);
    %     %T = aggregateROIsFromALF_fast(sPaths,'columnName',statName,'useChronic',true);
    %
    %     %generate density maps
    %     fig = IBL_plotROIdensity(T,statName);
    %
    %     %save table and figure
    %     if saveflag
    %         parquetwrite(fullfile(svPath,['allROIs_',statName]),T);
    %         saveas(fig,fullfile(svPath,['allROIs_',statName]),'fig');
    %         saveas(fig,fullfile(svPath,['allROIs_',statName]),'png');
    %     end
    %
    % end
end

if false
%% passive movie reliability map

% Aggregate once (any valid columnName will do; we only use passiveMovieCorr)
Tmovie = aggregateROIsFromALF_fast( ...
    sPaths, ...
    'useType', true);

D = IBL_computeROIdensityMaps(Tmovie);
D.v_clim_sess = [0 1];
%cmapPassive = brewermap(256,'Reds');
cmapPassive = parula;

nSess = height(unique(Tmovie(:,1:3),'rows'));
nSubj = height(unique(Tmovie(:,1),'rows'));

figMovie = figure( ...
    'Position', [1500 -200 500 350], ...
    'Name', 'PassiveMovieReliability', ...
    'Color', 'k');

tlMovie = tiledlayout(1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

axMovie = nexttile(tlMovie, 1);

IBL_plotStatPanel( ...
    D, ...
    D.passiveMovieCorr_mean, ...
    D.passiveMovieCorr_clim, ...
    'MapType','sequential',...
    'Colormap',cmapPassive);

title(sprintf('Passive movie reliability (n=%d)',nSess), ...
    'Color','w', ...
    'FontWeight','bold');

axLeg = nexttile(tlMovie, 2);

plotStatLegend( ...
    axLeg, ...
    D.passiveMovieCorr_clim, ...
    D.v_clim_sess, ...
    'mean repeat corr.', ...
    'MapType', 'sequential',...
    'Colormap',cmapPassive);

if saveflag
    exportgraphics( ...
        figMovie, ...
        fullfile(svPath,'ROIdensity_PassiveMovieReliability.png'), ...
        'Resolution',300, ...
        'BackgroundColor','current');
end

%% ROI summary map


%make custom colormap
cmap = brewermap(256,'Reds');
cmap(1,:) = [0 0 0];
nRamp = 80;
for i = 2:nRamp
    t = (i-1)/(nRamp-1);
    cmap(i,:) = t*cmap(nRamp,:);
end

%reduce table to one entry per cUID
Tmovie_red = reduceROITableByCUID(Tmovie);

%subsample to approach spatially uniform
Tmovie_red_unif = spatialSubsampleROIs(Tmovie,...
    'binSize',100,'maxPerBin',1000,...
    'verbose',true);

figROIMovie = figure( ...
    'Position', [1500 -200 700 700], ...
    'Name', 'ROIPassiveMovieCorrMap', ...
    'Color', 'k');

axROIMovie = axes(figROIMovie, 'Color', 'k');

IBL_plotROIPassiveMovieCorrMap( ...
    Tmovie_red_unif, ...
    'ax', axROIMovie, ...
    'clim', [0 .6], ...
    'dotSize', 4, ...
    'dotAlpha', 0.3, ...
    'showColorbar', true, ...
    'colormap',cmap);

title(axROIMovie, ...
    'Passive movie cross-repeat reliability', ...
    'Color', 'w', ...
    'Interpreter', 'none');

if saveflag
    exportgraphics( ...
        figROIMovie, ...
        fullfile(svPath,'ROISummaryMap_PassiveMovieReliability.png'), ...
        'Resolution',300, ...
        'BackgroundColor','current');
end

%% Passive stimulus responsive-ROI fraction maps: valve, tone, and noise

passiveStats = {'valve', 'tone', 'noise'};

passiveTitles = { ...
    'Valve-responsive ROIs', ...
    'Tone-responsive ROIs', ...
    'Noise-responsive ROIs'};

passiveLabels = { ...
    'fraction valve responsive', ...
    'fraction tone responsive', ...
    'fraction noise responsive'};

cmapPassiveStats = parula;

for iStat = 1:numel(passiveStats)

    statName = passiveStats{iStat};

    %Load ROI-level passive tuning results

    Tpassive = aggregateROIsFromALF_fast( ...
        sPaths, ...
        'columnName', statName, ...
        'useType', true, ...
        'alpha', 0.05);

    if isempty(Tpassive) || height(Tpassive) == 0
        warning( ...
            'No ROI data found for passive statistic "%s".', ...
            statName);
        continue
    end


    %Compute spatial maps

    D = IBL_computeROIdensityMaps(Tpassive);

    % Fraction of significant ROIs in each spatial bin.
    % For passive tuning, Tpassive.h is true when p < 0.05.
    responsiveFractionMap = D.ratio;

    % Plotting limits.
    D.frac_clim = [0 0.2];
    D.v_clim_sess = [0 1];


    %Session and subject counts

    nSess = height(unique( ...
        Tpassive(:, {'subject', 'date', 'session'}), ...
        'rows'));

    nSubj = height(unique( ...
        Tpassive(:, {'subject'}), ...
        'rows'));


    %Create figure

    figPassive = figure( ...
        'Position', [1500 -200 500 350], ...
        'Name', sprintf('Passive_%s_ResponsiveFraction', statName), ...
        'Color', 'k');

    tlPassive = tiledlayout( ...
        figPassive, ...
        1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');


    %Responsive-fraction map

    axPassive = nexttile(tlPassive, 1);

    IBL_plotStatPanel( ...
        D, ...
        responsiveFractionMap, ...
        D.frac_clim, ...
        'MapType', 'sequential', ...
        'Colormap', cmapPassiveStats);

    title( ...
        sprintf( ...
            '%s\n(n=%d sessions, %d subjects)', ...
            passiveTitles{iStat}, ...
            nSess, ...
            nSubj), ...
        'Color', 'w', ...
        'FontWeight', 'bold');


    %Legend

    axLeg = nexttile(tlPassive, 2);

    plotStatLegend( ...
        axLeg, ...
        D.frac_clim, ...
        D.v_clim_sess, ...
        passiveLabels{iStat}, ...
        'MapType', 'sequential', ...
        'Colormap', cmapPassiveStats);


    %Save figure

    if saveflag
        exportgraphics( ...
            figPassive, ...
            fullfile( ...
                svPath, ...
                sprintf( ...
                    'ROIdensity_Passive%sResponsiveFraction.png', ...
                    upperFirst(statName))), ...
            'Resolution', 300, ...
            'BackgroundColor', 'current');
    end

end

end

%% make session-by-session heatmaps per stat/base
if false
    for b = 1:numel(bases)

        bparts = bases{b};
        saveBase = sprintf('%s_%s_%s', bparts{1}, bparts{2}, bparts{3});
        outBase = ['ROIdensity_sessions_', saveBase];

        nTime = numel(twin_ev);

        % Store density maps for each time bin
        D_all = cell(nTime,1);
        sessKeys_all = cell(nTime,1);
        sessLabels_all = cell(nTime,1);

        for i = 1:nTime

            t = twin_ev{i} * 1000;
            t_mid = round(mean(t));

            statName = sprintf('%s_%s_%dto%d_%s', ...
                bparts{1}, bparts{2}, t(1), t(2), bparts{3});

            parquetPath = fullfile(svPath, ['allROIs_', statName, '.parquet']);

            if loadParquetIfExists && ~forceReaggregate && isfile(parquetPath)

                fprintf('Loading saved ROI table: %s\n', parquetPath);
                T = parquetread(parquetPath);

            else

                fprintf('Aggregating ROI table from ALF: %s\n', statName);
                T = aggregateROIsFromALF_fast( ...
                    sPaths, ...
                    'columnName', statName, ...
                    'useType', true, ...
                    'onlyResponsive', true);

                if saveflag
                    if ~exist(svPath, 'dir')
                        mkdir(svPath);
                    end

                    fprintf('Saving ROI table: %s\n', parquetPath);
                    parquetwrite(parquetPath, T);
                end

            end

            % Optional chronic reduction / spatial subsampling
            % T_red = reduceROITableByCUID(T);
            % T_red_unif = spatialSubsampleROIs(T_red);

            % Keep only cUIDs that are non-empty and present in every session
            % sessKeys = strcat(string(T.subject), "_", string(T.date), "_", string(T.session));
            % cUIDs = string(T.cUID);
            % hasCUID = cUIDs ~= "";
            % uSess = unique(sessKeys, 'stable');
            % nSess = numel(uSess);
            % cUIDs = string(T.cUID);
            % hasCUID = cUIDs ~= "";
            % Tc = table(cUIDs(hasCUID), sessKeys(hasCUID), ...
            %     'VariableNames', {'cUID','sessKey'});
            % [G, uidList] = findgroups(Tc.cUID);
            % nSessPerUID = splitapply(@(x) numel(unique(x)), Tc.sessKey, G);
            % sharedCUIDs = uidList(nSessPerUID == nSess);
            % T = T(hasCUID & ismember(cUIDs, sharedCUIDs), :);

            % Unique sessions in this (filtered) aggregate table
            sessKeys = strcat(string(T.subject), "_", string(T.date), "_", string(T.session));
            [uSessKeys, ia] = unique(sessKeys, 'stable');
            sessLabels = strcat( ...
                string(T.subject(ia)), " | ", ...
                string(T.date(ia)), " | ", ...
                string(T.session(ia)));
            nSess = numel(uSessKeys);

            D_sess = cell(nSess,1);

            for s = 1:nSess
                Ts = T(sessKeys == uSessKeys(s), :);

                D = IBL_computeROIdensityMaps(Ts);
                D.v_clim_sess = [0 1];

                D_sess{s} = D;
            end

            D_all{i} = D_sess;
            sessKeys_all{i} = uSessKeys;
            sessLabels_all{i} = sessLabels;

        end

        % Use union of sessions across all time bins
        allSessKeys = unique(vertcat(sessKeys_all{:}), 'stable');
        nSessTotal = numel(allSessKeys); %min([4,numel(allSessKeys)]);

        % ----- One color scale per session across all time bins and panel types -----
        sessCLim = containers.Map;
        for s = 1:nSessTotal
            vals = [];
            for i = 1:nTime
                hit = find(sessKeys_all{i} == allSessKeys(s), 1);
                if isempty(hit)
                    continue;
                end
                D = D_all{i}{hit};
                vals = [vals; D.ratio_hi(:); D.ratio_lo(:)]; %#ok<AGROW>
            end
            vals = vals(isfinite(vals));
            if isempty(vals)
                sessCLim(char(allSessKeys(s))) = [0 1];
            else
                lo = prctile(vals, 10);
                hi = prctile(vals, 99);
                if lo == hi
                    hi = lo + eps;
                end
                sessCLim(char(allSessKeys(s))) = [lo hi];
            end
        end

        allSessLabels = strings(nSessTotal,1);
        for s = 1:nSessTotal
            for i = 1:nTime
                hit = find(sessKeys_all{i} == allSessKeys(s), 1);
                if ~isempty(hit)
                    allSessLabels(s) = sessLabels_all{i}(hit);
                    break;
                end
            end
        end

        panelTypes = { ...
            'overall', 'ratio',    'tot frac. sig.'; ...
            'high',    'ratio_hi', 'hi frac. sig.'; ...
            'low',     'ratio_lo', 'lo frac. sig.'};

        for pIdx = 1:size(panelTypes,1)

            panelName  = panelTypes{pIdx,1};
            panelField = panelTypes{pIdx,2};
            panelLabel = panelTypes{pIdx,3};

            fig = figure( ...
                'Position',[1500,-300,1200,max(400,260*nSessTotal)], ...
                'Name',[outBase, '_', panelName], ...
                'Color','k');

            nCols = nTime + 1;

            tl = tiledlayout(nSessTotal, nCols, ...
                'TileSpacing','compact', ...
                'Padding','compact');

            for s = 1:nSessTotal

                for i = 1:nTime

                    t = twin_ev{i} * 1000;
                    t_mid = round(mean(t));
                    colLabel = sprintf('%d ms', t_mid);

                    ax = nexttile((s-1)*nCols + i);

                    hit = find(sessKeys_all{i} == allSessKeys(s), 1);

                    if ~isempty(hit)

                        D = D_all{i}{hit};
                        IBL_plotDensityPanel(D, D.(panelField), cmaps_frac_perBase{b});
                        %IBL_plotDensityPanel(D, D.(panelField), sessCLim(char(allSessKeys(s))));

                    else

                        axis(ax, 'off');
                        set(ax, 'Color', 'k');

                    end

                    if s == 1
                        title(ax, colLabel, 'Color','w');
                    end

                    if i == 1
                        ylabel(ax, allSessLabels(s), ...
                            'Color','w', ...
                            'FontWeight','bold', ...
                            'Interpreter','none');
                    end

                end

                % right-hand-side HSV legend for this session row
                axLeg = nexttile((s-1)*nCols + nCols);
                plotHSVLegend( ...
                    axLeg, ...
                    cmaps_frac_perBase{b},... %sessCLim(char(allSessKeys(s))), ...
                    D.v_clim_sess, ...
                    panelLabel);

                % % Empty last-column tile except for one legend
                % axLegSlot = nexttile((s-1)*nCols + nCols);
                % axis(axLegSlot, 'off');
                % set(axLegSlot, 'Color','k');

            end

            % % right-hand-side HSV legend
            % axLeg = nexttile(nCols);
            % plotHSVLegend( ...
            %     axLeg, ...
            %     cmaps_frac_perBase{b}, ...
            %     [0 1], ...
            %     panelLabel);

            sgtitle(sprintf('%s | %s', baseLabels{b}, panelLabel), ...
                'Color','w', ...
                'Interpreter','none');

            if saveflag
                if ~exist(svPath, 'dir')
                    mkdir(svPath);
                end

                exportgraphics(fig, ...
                    fullfile(svPath, sprintf('%s_%s.png', outBase, panelName)), ...
                    'Resolution', 300, ...
                    'BackgroundColor', 'current');
            end

        end

    end
end
%% make session-by-session ROI scatter figures per stat/base
if false %DRAFT
    for b = 1:numel(bases)

        bparts = bases{b};
        saveBase = sprintf('%s_%s_%s', bparts{1}, bparts{2}, bparts{3});
        outBase = ['ROISummaryMap_sessions_', saveBase];

        nTime = numel(twin_ev);

        T_all = cell(nTime,1);
        sessKeys_all = cell(nTime,1);
        sessLabels_all = cell(nTime,1);

        % ----- Load / aggregate all time bins -----
        for i = 1:nTime

            t = twin_ev{i} * 1000;
            t_mid = round(mean(t));

            statName = sprintf('%s_%s_%dto%d_%s', ...
                bparts{1}, bparts{2}, t(1), t(2), bparts{3});

            parquetPath = fullfile(svPath, ['allROIs_', statName, '.parquet']);

            if loadParquetIfExists && ~forceReaggregate && isfile(parquetPath)

                fprintf('Loading saved ROI table: %s\n', parquetPath);
                T = parquetread(parquetPath);

            else

                fprintf('Aggregating ROI table from ALF: %s\n', statName);
                T = aggregateROIsFromALF_fast( ...
                    sPaths, ...
                    'columnName', statName, ...
                    'useType', true, ...
                    'onlyResponsive', true, ...
                    'useChronic', true);

                if saveflag
                    if ~exist(svPath, 'dir')
                        mkdir(svPath);
                    end

                    fprintf('Saving ROI table: %s\n', parquetPath);
                    parquetwrite(parquetPath, T);
                end

            end

            T.cUID = string(T.cUID);

            sessKeys = strcat(string(T.subject), "_", string(T.date), "_", string(T.session));

            T_all{i} = T;
            sessKeys_all{i} = sessKeys;

            [uSessKeys, ia] = unique(sessKeys, 'stable');
            sessLabels_all{i} = table( ...
                uSessKeys, ...
                strcat(string(T.subject(ia)), " | ", string(T.date(ia)), " | ", string(T.session(ia))), ...
                'VariableNames', {'sessKey','label'});

        end

        % ----- Shared sessions across time bins -----
        allSessKeys = sessKeys_all{1};
        allSessKeys = unique(allSessKeys, 'stable');

        for i = 2:nTime
            allSessKeys = allSessKeys(ismember(allSessKeys, unique(sessKeys_all{i}, 'stable')));
        end

        nSessTotal = numel(allSessKeys);

        allSessLabels = strings(nSessTotal,1);
        for s = 1:nSessTotal
            for i = 1:nTime
                L = sessLabels_all{i};
                hit = find(L.sessKey == allSessKeys(s), 1);
                if ~isempty(hit)
                    allSessLabels(s) = L.label(hit);
                    break;
                end
            end
        end

        % ----- Keep only cUIDs shared across all sessions and all time bins -----
        sharedCUIDs = [];

        for i = 1:nTime

            T = T_all{i};
            sessKeys = sessKeys_all{i};

            hasCUID = T.cUID ~= "";

            perTimeShared = string.empty(0,1);

            for s = 1:nSessTotal
                thisUIDs = unique(T.cUID(hasCUID & sessKeys == allSessKeys(s)));
                if s == 1
                    perTimeShared = thisUIDs;
                else
                    perTimeShared = perTimeShared(ismember(perTimeShared, thisUIDs));
                end
            end

            if i == 1
                sharedCUIDs = perTimeShared;
            else
                sharedCUIDs = sharedCUIDs(ismember(sharedCUIDs, perTimeShared));
            end

        end

        fprintf('%s: keeping %d cUIDs shared across all sessions/time bins.\n', ...
            saveBase, numel(sharedCUIDs));

        % ----- Reduce + spatially subsample each time/session table -----
        T_plot = cell(nTime, nSessTotal);

        for i = 1:nTime

            T = T_all{i};
            sessKeys = sessKeys_all{i};

            keepShared = T.cUID ~= "" & ismember(T.cUID, sharedCUIDs);
            T = T(keepShared,:);

            for s = 1:nSessTotal

                Ts = T(sessKeys(keepShared) == allSessKeys(s), :);

                if isempty(Ts)
                    T_plot{i,s} = Ts;
                    continue;
                end

                % Reduce to one row per cUID within this session/time bin
                Ts_red = reduceROITableByCUID(Ts);

                % Spatially subsample to reduce overplotting / flatten density
                Ts_red_unif = spatialSubsampleROIs(Ts_red,'verbose',false);

                T_plot{i,s} = Ts_red_unif;

            end
        end

        % ----- Auto limits from full aggregate tables -----
        allML = [];
        allAP = [];

        for i = 1:nTime
            allML = [allML; T_all{i}.ML]; %#ok<AGROW>
            allAP = [allAP; T_all{i}.AP]; %#ok<AGROW>
        end

        validPos = isfinite(allML) & isfinite(allAP);

        pad = 100; % microns

        xLimits = [min(allML(validPos)) max(allML(validPos))] + [-pad pad];
        yLimits = [min(allAP(validPos)) max(allAP(validPos))] + [-pad pad];

        % Optional: round to nearest 50 um
        xLimits = 50 * [floor(xLimits(1)/50), ceil(xLimits(2)/50)];
        yLimits = 50 * [floor(yLimits(1)/50), ceil(yLimits(2)/50)];

        % ----- Plot one figure: columns = time bins, rows = sessions -----
        figMap = figure( ...
            ...%'Position',[1500,-200,800,max(300,180*nSessTotal)], ...
            'Position',[1500,-200,1500,max(500,280*nSessTotal)], ...
            'Name',outBase, ...
            'Color','k');

        tl = tiledlayout(nSessTotal, nTime, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        for s = 1:nSessTotal

            for i = 1:nTime

                t = twin_ev{i} * 1000;
                t_mid = round(mean(t));
                colLabel = sprintf('%d ms', t_mid);

                ax = nexttile((s-1)*nTime + i);

                Ts = T_plot{i,s};

                if ~isempty(Ts)

                    IBL_plotROISummaryMap( ...
                        Ts, ...
                        'ax', ax, ...
                        'alpha2', 0.025, ...
                        'plotNonsigAsDots', true,...
                        'xLimits', xLimits, ...
                        'yLimits', yLimits);

                else

                    axis(ax, 'off');
                    set(ax, 'Color','k');

                end

                if s == 1
                    title(ax, colLabel, 'Color','w');
                end

                if i == 1
                    ylabel(ax, allSessLabels(s), ...
                        'Color','w', ...
                        'FontWeight','bold', ...
                        'Interpreter','none');
                end

            end

        end

        sgtitle(sprintf('%s | shared cUID ROI summary maps', baseLabels{b}), ...
            'Color','w', ...
            'Interpreter','none');

        if saveflag
            if ~exist(svPath, 'dir')
                mkdir(svPath);
            end

            exportgraphics(figMap, ...
                fullfile(svPath, [outBase '.png']), ...
                'Resolution', 300, ...
                'BackgroundColor', 'current');
        end

        if b==1

        figLeg = figure( ...
            'Position',[1500,-200,170,170], ...
            'Name',['ROISummaryLegend_', saveBase], ...
            'Color','k');

        axLeg = axes(figLeg, 'Color','k');

        IBL_plotROISummaryLegend( ...
            axLeg, ...
            'minSize', 2, ...
            'sizeScale', 0.05, ...
            'alpha2', 0.025,...
            'label',bparts{1});

        if saveflag
            exportgraphics(figLeg, ...
                fullfile(svPath, ['ROISummaryLegend_', saveBase, '.png']), ...
                'Resolution', 300, ...
                'BackgroundColor', 'current');
        end

        figInset = figure( ...
            'Position',[1500,-200,250,300], ...
            'Name',['ROISummaryMap_inset_', saveBase], ...
            'Color','k');

        axInset = axes(figInset, 'Color','k');
        hold(axInset,'on');

        plotROILimitInset(axInset, xLimits, yLimits);

        if saveflag
            exportgraphics(figInset, ...
                fullfile(svPath, ['ROISummaryMap_inset_', saveBase, '.png']), ...
                'Resolution', 300, ...
                'BackgroundColor', 'current');
        end

        end

    end
end


%% helpers
function plotHSVLegend(ax, fraclim, vlim, labelstr)

if nargin < 4
    labelstr = 'frac. sig.';
end

n = 250;

[V, F] = meshgrid(linspace(vlim(1),   vlim(2),   n), ...
    linspace(fraclim(1), fraclim(2), n));

Fnorm = min(max((F - fraclim(1)) ./ diff(fraclim), 0), 1);
Vnorm = min(max((V - vlim(1)) ./ diff(vlim), 0), 1);

H = 0.70 * (1 - Fnorm);
S = ones(size(H));

RGB = hsv2rgb(cat(3, H, S, Vnorm));

image(ax, ...
    linspace(vlim(1),   vlim(2),   n), ...
    linspace(fraclim(1), fraclim(2), n), ...
    RGB);

set(ax, ...
    'YDir','normal', ...
    'XColor','w', ...
    'YColor','w', ...
    'Color','k', ...
    'Box','on');

xlabel(ax, 'nr. sess', 'FontSize', 10, 'Color',[1 1 1]);
ylabel(ax, labelstr,  'FontSize', 10, 'Color',[1 1 1]);

axis(ax,'tight');
%axis(ax,'square');
ax.PlotBoxAspectRatio = [1,5,1];

end

function plotStatLegend(ax, clim, vlim, labelstr, varargin)
%PLOTSTATLEGEND Plot a 2-D color/brightness legend.
%
% Color encodes the statistic value.
% Brightness encodes the number of contributing sessions.
%
% Name-value options:
%   'MapType'  : 'diverging' or 'sequential'
%   'Center'   : center of diverging map, default 0
%   'Colormap' : explicit n-by-3 colormap

if nargin < 4 || isempty(labelstr)
    labelstr = 'stat';
end

p = inputParser;
p.addParameter('MapType', 'diverging', ...
    @(x) ischar(x) || isstring(x));
p.addParameter('Center', 0, ...
    @(x) isnumeric(x) && isscalar(x));
p.addParameter('Colormap', [], ...
    @(x) isempty(x) || (isnumeric(x) && size(x,2) == 3));
p.parse(varargin{:});

mapType = validatestring( ...
    char(p.Results.MapType), ...
    {'diverging','sequential'});

center = p.Results.Center;

n = 250;
nC = 256;

% Resolve color limits
switch mapType
    case 'diverging'
        if numel(clim) == 1
            d = abs(clim);
        elseif numel(clim) == 2
            d = max(abs(clim - center));
        else
            error('Diverging clim must contain one or two values.');
        end

        clim = center + [-d d];

    case 'sequential'
        if numel(clim) == 1
            clim = [0 clim];
        elseif numel(clim) ~= 2
            error('Sequential clim must contain one or two values.');
        end
end

if clim(2) <= clim(1)
    error('clim must be increasing.');
end

if vlim(2) <= vlim(1)
    error('vlim must be increasing.');
end

% Choose color map
if isempty(p.Results.Colormap)
    switch mapType
        case 'diverging'
            cmap = brewermap(nC, '*RdBu');
        case 'sequential'
            cmap = brewermap(nC, 'YlGnBu');
    end
else
    cmap = p.Results.Colormap;
    nC = size(cmap,1);
end

% Grid: horizontal = session count, vertical = statistic
[V, S] = meshgrid( ...
    linspace(vlim(1), vlim(2), n), ...
    linspace(clim(1), clim(2), n));

% Statistic controls color
Snorm = (S - clim(1)) ./ diff(clim);
Snorm = min(max(Snorm, 0), 1);

idx = round(1 + Snorm .* (nC - 1));
RGB = ind2rgb(idx, cmap);

% Session count controls brightness
Vnorm = (V - vlim(1)) ./ diff(vlim);
Vnorm = min(max(Vnorm, 0), 1);

RGB = RGB .* Vnorm;

% Plot
image( ...
    ax, ...
    linspace(vlim(1), vlim(2), n), ...
    linspace(clim(1), clim(2), n), ...
    RGB);

set(ax, ...
    'YDir', 'normal', ...
    'XColor', 'w', ...
    'YColor', 'w', ...
    'Color', 'k', ...
    'Box', 'on');

xlabel(ax, 'nr. sess', ...
    'FontSize', 10, ...
    'Color', 'w');

ylabel(ax, labelstr, ...
    'FontSize', 10, ...
    'Color', 'w', ...
    'Interpreter', 'tex');

axis(ax, 'tight');
ax.PlotBoxAspectRatio = [1 5 1];

end

function IBL_plotROISummaryLegend(ax, varargin)

p = inputParser;
addParameter(p, 'minSize', 1, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'sizeScale', 0.02, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'statExamples', [10 500 2000], @(x)isnumeric(x));
addParameter(p, 'alpha2', 0.025, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'label', 'stat', @(x)ischar(x) || isstring(x));
parse(p, varargin{:});
opts = p.Results;

cla(ax);
hold(ax,'on');
axis(ax,[0 1 0 1]);
axis(ax,'off');
set(ax,'Color','k');

cBlue = [0.35 0.60 0.95];
cRed  = [0.95 0.45 0.45];
cGray = [0.55 0.55 0.55];

colors = [
    cBlue
    cGray
    cRed
    ];

pLabels = {
    sprintf('p < %.3g', opts.alpha2)
    'n.s.'
    sprintf('p > %.3g', 1 - opts.alpha2)
    };

statExamples = opts.statExamples(:)';
nStat = numel(statExamples);

if nStat ~= 3
    error('statExamples should contain exactly 3 values for this 3 x 3 legend.');
end

x = [0.28 0.52 0.76];
y = [0.72 0.50 0.28];

% title / axis labels
% text(ax,0.52,0.94,'ROI dot legend', ...
%     'Color','w', ...
%     'FontWeight','bold', ...
%     'HorizontalAlignment','center');

text(ax,0.52,0.06,sprintf('|%s|', opts.label), ...
    'Color','w', ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center');

% text(ax,0.07,0.50,'pval', ...
%     'Color','w', ...
%     'FontWeight','bold', ...
%     'HorizontalAlignment','center', ...
%     'Rotation',90);

% dots
for r = 1:3
    for c = 1:3
        st = statExamples(c);
        sz = opts.minSize + opts.sizeScale * abs(st);

        scatter(ax, x(c), y(r), sz, colors(r,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.8, ...
            'MarkerEdgeAlpha', 0.2);
    end
end

% x tick labels: stat examples
for c = 1:3
    text(ax, x(c), 0.22, sprintf('%g', statExamples(c)), ...
        'Color','w', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top',...
        'Fontsize',8);
end

% y tick labels: p-value categories
for r = 1:3
    text(ax, 0.2, y(r), pLabels{r}, ...
        'Color','w', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle',...
        'Fontsize',8);
end

end

function plotROILimitInset(ax, xLimits, yLimits, varargin)

p = inputParser;
addParameter(p, 'bas', aratopdown.atlas.build_topdown, @isstruct);
parse(p, varargin{:});
opts = p.Results;

cla(ax);
hold(ax,'on');

set(ax, ...
    'Color','k', ...
    'XColor','none', ...
    'YColor','none', ...
    'XTick',[], ...
    'YTick',[], ...
    'Box','off');

% Full atlas boundaries
cellfun(@(x) cellfun(@(y) ...
    plot(ax, 1000*y(:,2), 1000*y(:,1), ...
    'Color',[0.45 0.45 0.45], ...
    'LineWidth',0.5), ...
    x, 'uni', false), ...
    {opts.bas.dorsal_brain_areas(1:end-11).boundaries_stereotax}, ...
    'uni', false);

% ROI map window
rectangle(ax, ...
    'Position', [xLimits(1), yLimits(1), diff(xLimits), diff(yLimits)], ...
    'EdgeColor', 'w', ...
    'LineWidth', 1.5);

axis(ax,'equal');
axis(ax,'tight');

end


function outputText = upperFirst(inputText)
%UPPERFIRST Capitalize the first character of a character vector.

inputText = char(inputText);

if isempty(inputText)
    outputText = inputText;
else
    outputText = [upper(inputText(1)), inputText(2:end)];
end

end
