
%datpath = 'Y:\Subjects\SP044\2023-06-27\001'; %bad one
%datpath = 'Y:\Subjects\SP044\2023-06-30\001';
%datpath = 'Y:\Subjects\SP044\2023-07-03\001';
%datpath = 'Y:\Subjects\SP081\2026-09-08\001'; %good one
%datpath = 'Y:\Subjects\SP081\2026-09-10\001'; %illustrate ROI saturation?
datpath = 'Y:\Subjects\SP076\2025-11-10\001'; %bad?

fn_trace = 'mpci.ROIActivityDeconvolved.npy';
fn_F = 'mpci.ROIActivityF.npy';
fn_Fneu = 'mpci.ROINeuropilActivityF.npy';
fn_meanImg = 'mpciMeanImage.images.npy';
fn_times = 'mpci.times.npy';
fn_badframes = 'mpci.badFrames.npy';
fn_iscell = 'mpciROIs.cellClassifier.npy';

classThresh = 0.5;


%% load and plot mean images

ff = fullfile(datpath,'alf');

paths_meta = dir(fullfile(datpath,'**','*000*.mat'));
load(fullfile(paths_meta(1).folder,paths_meta(1).name));

paths_FOV = dir(fullfile(datpath,'alf','FOV_*',fn_trace));
nFOVs = length(paths_FOV);

meanImgs = {};
for ii=1:nFOVs
    meanImgs{ii} = readNPY([paths_FOV(ii).folder filesep fn_meanImg]);
end

% Common robust color limits
allPixels = cell2mat(cellfun(@(x) x(:), meanImgs, 'UniformOutput', false));
colorLimits = prctile(allPixels(:), [1 99]);

figure;
for ii=1:nFOVs
    subplot(2,ceil(nFOVs/2),ii);
    imagesc(meanImgs{ii});
    clim(colorLimits);
    colormap gray;
    title(sprintf('FOV_%02d, z=%d',ii-1,meta.FOV(ii).Zs),'Interpreter','none');
    axis image;
    axis off;
end
sgtitle(datpath,'Interpreter','none');

%% compare traces from two FOVs

FOVs_toplot = {'FOV_01','FOV_05'};
tr_all = {};

for i=1:length(FOVs_toplot)
    tr_all{i} =  readNPY([ff filesep FOVs_toplot{i} filesep fn_trace]);
end

figure;
for i=1:length(FOVs_toplot)
    subplot(length(FOVs_toplot),1,i);
    plot(tr_all{i}(1:1000,1));
    title(sprintf('%s, ROI 0',FOVs_toplot{i}),'Interpreter','none');
    xlabel('frame')
    ylabel(fn_trace(1:end-4));
end
sgtitle(datpath,'Interpreter','none');

%% compute neural quality metrics

qmNames = {'noiseLevel','mean','std','skew','saturationRatio'};

neuralQMs_all = cell(1,nFOVs);
iscell_all = cell(1,nFOVs);

for ii = 1:nFOVs

    fovDir = paths_FOV(ii).folder;

    F = readNPY(fullfile(fovDir,fn_F));
    Fneu = readNPY(fullfile(fovDir,fn_Fneu));
    times = readNPY(fullfile(fovDir,fn_times));

    % Bad frames
    badframesPath = fullfile(fovDir,fn_badframes);
    if exist(badframesPath,'file')
        badframes = logical(readNPY(badframesPath));
    else
        badframes = false(size(F,1),1);
    end
    badframes = badframes(:);

    % Cell classifier
    classPath = fullfile(fovDir,fn_iscell);
    if exist(classPath,'file')
        classVals = readNPY(classPath);
        iscell = classVals>classThresh;
        if size(iscell,2) > 1
            iscell = iscell(:,1);
        end
    else
        iscell = true(size(F,2),1);
    end
    iscell = iscell(:);

    iscell_all{ii} = iscell;

    neuralQMs_all{ii} = get_neuralQMs( ...
        F,Fneu,times, ...
        'badframes',badframes, ...
        'iscell',iscell);

end


%% mean images + neural quality metric distributions

nQM = numel(qmNames);

figure;
t = tiledlayout(nQM+1,nFOVs, ...
    'TileSpacing','loose', ...
    'Padding','compact');

%% top row: mean images

for ii = 1:nFOVs

    ax = nexttile(ii);

    imagesc(meanImgs{ii});
    clim(colorLimits);
    colormap(ax,gray);

    axis image off;

    title(sprintf('FOV_%02d\nz=%d',ii-1,meta.FOV(ii).Zs), ...
        'Interpreter','none');

    nCell = sum(iscell_all{ii});
    nNonCell = sum(~iscell_all{ii});

    text(ax,0.5,-0.04, ...
        sprintf('cell: %d   non-cell: %d',nCell,nNonCell), ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',8);

end


%% quality metrics

for iQM = 1:nQM

    qmName = qmNames{iQM};

    % Pool values across all FOVs to define common vertical range
    valsAll = [];

    for ii = 1:nFOVs
        vals = [neuralQMs_all{ii}.(qmName)];
        valsAll = [valsAll vals(isfinite(vals))];
    end

    ylims = prctile(valsAll,[1 99]);

    % Small padding so distributions / labels do not touch limits
    dy = diff(ylims);
    ylims = ylims + [-0.05 0.12]*dy;

    for ii = 1:nFOVs

        tileIdx = iQM*nFOVs + ii;
        ax = nexttile(tileIdx);
        hold(ax,'on');

        vals = [neuralQMs_all{ii}.(qmName)]';
        iscell = iscell_all{ii};

        valsCell = vals(iscell & isfinite(vals));
        valsNonCell = vals(~iscell & isfinite(vals));

        nCell = numel(valsCell);
        nNonCell = numel(valsNonCell);

        % Common evaluation points for KDE
        yEval = linspace(ylims(1),ylims(2),150);

        % Non-cells: LEFT, gray
        if nNonCell > 1
            fNonCell = ksdensity(valsNonCell,yEval);
            fNonCell = fNonCell ./ max(fNonCell);

            plot(ax,-fNonCell,yEval, ...
                'Color',[0.6 0.6 0.6], ...
                'LineWidth',1.5);
        end

        % Cells: RIGHT, colored
        if nCell > 1
            fCell = ksdensity(valsCell,yEval);
            fCell = fCell ./ max(fCell);

            plot(ax,fCell,yEval, ...
                'LineWidth',1.8);
        end

        % Central dividing line
        xline(ax,0,':','Color',[0.7 0.7 0.7]);

        % Medians
        medCell = median(valsCell,'omitnan');
        medNonCell = median(valsNonCell,'omitnan');

        if isfinite(medNonCell)
            plot(ax,[-1 0], ...
                [medNonCell medNonCell], ...
                '-', ...
                'Color',[0.6 0.6 0.6], ...
                'LineWidth',1);
        end

        if isfinite(medCell)
            plot(ax,[0 1], ...
                [medCell medCell], ...
                '-', ...
                'LineWidth',1);
        end

        % Median values at top
        text(ax,0.25,1.02, ...
            sprintf('%.3g',medNonCell), ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',7, ...
            'Color',[0.5 0.5 0.5]);

        text(ax,0.75,1.02, ...
            sprintf('%.3g',medCell), ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',7);


        ylim(ax,ylims);
        xlim(ax,[-1.05 1.05]);

        set(ax,'XTick',[]);

        if ii == 1
            ylabel(ax,qmName,'Interpreter','none');
        else
            set(ax,'YTickLabel',[]);
        end

        box(ax,'off');

    end
end

%title(t,datpath,'Interpreter','none');

%% example from an FOV

fov_toplot = 'FOV_01';

fovDir = fullfile(datpath,'alf',fov_toplot);
F = readNPY(fullfile(fovDir,fn_F));

metricVals = [neuralQMs_all{2}.saturationRatio];

plotQMquantileTraces(F,metricVals, ...
    'mask',iscell_all{2}, ...
    'metricName','saturationRatio', ...
    'traceName','raw fluorescence');

%% example from bad FOV
fov_toplot = 'FOV_03';

fovDir = fullfile(datpath,'alf',fov_toplot);
F = readNPY(fullfile(fovDir,fn_F));

metricVals = [neuralQMs_all{4}.saturationRatio];

plotQMquantileTraces(F,metricVals, ...
    'mask',iscell_all{7}, ...
    'metricName','saturationRatio', ...
    'traceName','raw fluorescence');

%% old
%%% compute and plot neural quality metrics
%
% qmNames = {'noiseLevel','mean','std','skew'};
%
% neuralQMs_all = cell(1,nFOVs);
% fovQMs_all = cell(1,nFOVs);
% iscell_all = cell(1,nFOVs);
%
% for ii = 1:nFOVs
%
%     fovDir = paths_FOV(ii).folder;
%
%     % Load raw fluorescence + neuropil
%     F = readNPY(fullfile(fovDir,fn_F));
%     Fneu = readNPY(fullfile(fovDir,fn_Fneu));
%
%     % Frame times
%     times = readNPY(fullfile(fovDir,fn_times));
%
%     % Bad frames
%     badframesPath = fullfile(fovDir,fn_badframes);
%     if exist(badframesPath,'file')
%         badframes = logical(readNPY(badframesPath));
%     else
%         badframes = false(size(F,1),1);
%     end
%     badframes = badframes(:);
%
%     % Cell classifier
%     iscellPath = fullfile(fovDir,fn_iscell);
%     if exist(iscellPath,'file')
%         iscell = logical(readNPY(iscellPath));
%         iscell = iscell(:);
%     else
%         iscell = true(size(F,2),1);
%     end
%
%     iscell_all{ii} = iscell;
%
%     [neuralQMs_all{ii}, fovQMs_all{ii}] = get_neuralQMs( ...
%         F,Fneu,times, ...
%         'badframes',badframes, ...
%         'iscell',iscell);
%
% end
%
%
% %% neural QM distributions by FOV
%
% % Number of histogram bins
% nBins = 30;
%
% figure;
%
% for iQM = 1:numel(qmNames)
%
%     qmName = qmNames{iQM};
%
%     % Pool CELL values across FOVs to determine common robust limits/bins
%     valsAll = [];
%
%     for ii = 1:nFOVs
%         vals = [neuralQMs_all{ii}.(qmName)];
%         vals = vals(iscell_all{ii});
%         valsAll = [valsAll vals(isfinite(vals))];
%     end
%
%     % Robust plotting range prevents a few outliers from ruining the plots
%     xlims = prctile(valsAll,[1 99]);
%     edges = linspace(xlims(1),xlims(2),nBins+1);
%
%     for ii = 1:nFOVs
%
%         subplot(numel(qmNames),nFOVs,(iQM-1)*nFOVs + ii);
%
%         vals = [neuralQMs_all{ii}.(qmName)];
%         vals = vals(iscell_all{ii});
%
%         histogram(vals,edges,'Normalization','probability');
%
%         xlim(xlims);
%
%         if ii == 1
%             ylabel(qmName,'Interpreter','none');
%         end
%
%         if iQM == 1
%             title(sprintf('FOV_%02d',ii-1),'Interpreter','none');
%         end
%
%         if iQM < numel(qmNames)
%             set(gca,'XTickLabel',[]);
%         end
%
%         box off;
%     end
% end
%
% sgtitle(sprintf('%s - neural quality metrics',datpath),'Interpreter','none');
%
%
% %% FOV-level quality metrics
%
% fovQMmat = nan(nFOVs,numel(qmNames));
%
% for ii = 1:nFOVs
%     for iQM = 1:numel(qmNames)
%         fovQMmat(ii,iQM) = fovQMs_all{ii}.(qmNames{iQM});
%     end
% end
%
% figure;
%
% for iQM = 1:numel(qmNames)
%
%     subplot(2,2,iQM);
%
%     plot(0:nFOVs-1,fovQMmat(:,iQM),'o-');
%
%     xlabel('FOV');
%     ylabel(qmNames{iQM},'Interpreter','none');
%
%     xticks(0:nFOVs-1);
%     xlim([-0.5 nFOVs-0.5]);
%
%     box off;
%
% end
%
% sgtitle(sprintf('%s - FOV quality metrics',datpath),'Interpreter','none');
%
%
% %% overlaid neural QM distributions
%
% nBins = 30;
%
% figure;
%
% for iQM = 1:numel(qmNames)
%
%     qmName = qmNames{iQM};
%
%     valsAll = [];
%     for ii = 1:nFOVs
%         vals = [neuralQMs_all{ii}.(qmName)];
%         vals = vals(iscell_all{ii});
%         valsAll = [valsAll vals(isfinite(vals))];
%     end
%
%     xlims = prctile(valsAll,[1 99]);
%     edges = linspace(xlims(1),xlims(2),nBins+1);
%
%     subplot(2,2,iQM);
%     hold on;
%
%     for ii = 1:nFOVs
%
%         vals = [neuralQMs_all{ii}.(qmName)];
%         vals = vals(iscell_all{ii});
%
%         histogram(vals,edges, ...
%             'Normalization','probability', ...
%             'DisplayStyle','stairs', ...
%             'LineWidth',1.2);
%
%     end
%
%     xlim(xlims);
%     xlabel(qmName,'Interpreter','none');
%     ylabel('fraction of ROIs');
%     title(qmName,'Interpreter','none');
%     box off;
%
% end
%
% legend(arrayfun(@(x)sprintf('FOV_%02d',x),0:nFOVs-1, ...
%     'UniformOutput',false), ...
%     'Location','bestoutside');
%
% sgtitle(sprintf('%s - neural quality metrics',datpath),'Interpreter','none');