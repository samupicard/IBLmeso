%sPaths = {'Y:\Subjects\SP072\2025-08-26\001'};
%FOV_sel = [3,4,5,6];

% sPaths = {'Y:\Subjects\SP072\2025-08-26\001',...
%     ...%'Y:\Subjects\SP072\2025-08-27\001',...
%     'Y:\Subjects\SP072\2025-08-28\001',...
%     'Y:\Subjects\SP072\2025-09-02\002',...
%     'Y:\Subjects\SP072\2025-09-03\001',...
%     'Y:\Subjects\SP072\2025-09-04\001',...
%     'Y:\Subjects\SP072\2025-09-05\001'...
%     };
% FOV_sel = [3,4,5,6];

sPaths = {'Y:\Subjects\SS010\2026-10-12\003'};
%FOV_sel = [2,6];
FOV_sel = [];

results = compareRepeatCorr(sPaths, 'FOVs', FOV_sel, 'ClassifierThreshold', 0.5);

%% plot map of reliability
figure;
scatter(results.s2p.pos(:,1),results.s2p.pos(:,2),10,results.s2p.r,'filled');
colormap(flipud(gray));
clim([0.1,0.4]);

%% explore overall activity

for iSess = 1

X_mask = [results(iSess).masknmf.traces.repeat1;
          results(iSess).masknmf.traces.repeat2];
quantile(X_mask(:),[0.5,0.99])

X_s2p  = [results(iSess).s2p.traces.repeat1;
          results(iSess).s2p.traces.repeat2];

%do SVD
X_mask = X_mask - mean(X_mask,1);
X_s2p  = X_s2p  - mean(X_s2p,1);

[U_mask,S_mask,V_mask] = svd(X_mask,'econ');
[U_s2p, S_s2p, V_s2p ] = svd(X_s2p ,'econ');

v_mask = diag(S_mask).^2;
v_s2p  = diag(S_s2p ).^2;

v_mask = v_mask/sum(v_mask);
v_s2p  = v_s2p/sum(v_s2p);

%get timecourse of first principal component
pc1_mask = U_mask(:,1) * S_mask(1,1);
pc1_s2p = U_s2p(:,1) * S_s2p(1,1);

%sort ROIs by loading on PC1
[~,iSort_mask] = sort(V_mask(:,1));
[~,iSort_s2p] = sort(V_s2p(:,1));

cols = colororder;

%plot activity, sorted by PC1
figure; 
tiledlayout(20,1); 
nexttile([19 1]); 
imagesc(T,1:size(X_s2p,2),X_s2p(:,iSort_s2p)'); 
colormap(flipud(gray)); 
ylabel('suite2p ROIs, sorted by PC1'); 
set(gca,'Xtick',[]); 
box off; 
nexttile; 
plot(T,pc1_s2p,'LineWidth',2,'Color',cols(1,:)); 
box off; 
xlabel('time(s)');
sgtitle(results(iSess).sPath);

figure; 
tiledlayout(20,1); 
nexttile([19 1]); 
imagesc(T,1:size(X_mask,2),X_mask(:,iSort_s2p)'); 
colormap(flipud(gray)); 
ylabel('maskNMF ROIs, sorted by PC1'); 
set(gca,'Xtick',[]); 
box off; 
nexttile; 
plot(T,pc1_mask,'LineWidth',2,'Color',cols(1,:)); 
box off; 
xlabel('time(s)');
sgtitle(results(iSess).sPath);


%plot eigenspectra
figure;
loglog(v_s2p,'LineWidth',1.5);
hold on;
loglog(v_mask,'LineWidth',1.5);
xlim([1,3000]);
ylim([1e-6,1e-1]);

xlabel('PC');
ylabel('Variance');
legend({'suite2p','maskNMF'});
title(results(iSess).sPath);
box off;
end

%% choose session and get matches 

k=20;
iSess = 1;

%find putative matches between maskNMF and suite2p traces
X = [results(iSess).s2p.raw.repeat1; results(iSess).s2p.raw.repeat2];
Y = [results(iSess).masknmf.raw.repeat1; results(iSess).masknmf.raw.repeat2];

[bestS2pIdx,bestCorr] = findBestColumnMatches(X,Y,...
    'XFov', results(iSess).s2p.fovs, ...
    'YFov', results(iSess).masknmf.fovs, ...
    'XPos', results(iSess).s2p.pos, ...
    'YPos', results(iSess).masknmf.pos, ...
    'MaxDist', 10);

%maskNMF
[~,iROI_maskNMF]=maxk(results(iSess).masknmf.r,k);
iROI_s2p = bestS2pIdx(iROI_maskNMF);

% only consider matches with a sufficiently strong trace-match correlation
iMatches = ~isnan(bestS2pIdx(:)) & bestCorr(:) > 0.95;

bestS2pIdx_matches = bestS2pIdx(iMatches);
matchCorr_matches  = bestCorr(iMatches);


%% plot traces and their distributions, for both methods in matched neurons

t = [results(iSess).time; results(iSess).time(2)+results(iSess).time(end)+results(iSess).time];

X_s2p = [results(iSess).s2p.traces.repeat1(:,bestS2pIdx_matches); results(iSess).s2p.traces.repeat1(:,bestS2pIdx_matches)]; 
X_mask = [results(iSess).masknmf.traces.repeat1(:,iMatches); results(iSess).masknmf.traces.repeat1(:,iMatches)]; 

F_s2p = [results(iSess).s2p.raw.repeat1(:,bestS2pIdx_matches); results(iSess).s2p.raw.repeat1(:,bestS2pIdx_matches)]; 
F_mask = [results(iSess).masknmf.raw.repeat1(:,iMatches); results(iSess).masknmf.raw.repeat1(:,iMatches)]; 

for i = 1:10 %choose a ROI index here

x_s2p = X_s2p(:,i); 
x_mask = X_mask(:,i);
f_s2p = F_s2p(:,i); 
f_mask = F_mask(:,i);

figure; 
tl = tiledlayout(2,4); 

axTrace(1) = nexttile(tl, [1,3]); 
hold on;
plot(t,f_s2p,'Color',[.7 .7 .7]);
plot(t,x_s2p,'Color',cols(1,:));
box off

nexttile(tl, [1,1]); hold on
plotZeroInflatedGammaFit(x_s2p, 50, cols(1,:),'UpperQuantile',.99);
box off

axTrace(2) = nexttile(tl, [1,3]); 
hold on;
plot(t,f_mask,'Color',[.7 .7 .7]);
plot(t,x_mask,'Color',cols(2,:));
box off
xlabel('time (s)');

nexttile(tl, [1,1]); hold on
plotZeroInflatedGammaFit(x_mask, 50, cols(2,:),'UpperQuantile',.99);
box off
xlabel('value');

linkaxes(axTrace,'x');
xlim(axTrace,[0,200]);

axTrace(1).XAxis.Color = 'none';
axTrace(1).XTick = [];

end


%% compare ZIG parameters

zig_s2p  = fitZIGMatrix(X_s2p);
zig_mask = fitZIGMatrix(X_mask);

params = {'q','k','a','loc'};
labels = {'q = P(nonzero)', 'k shape', 'a scale', 'loc'};

figure;
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for iParam = 1:2

    p = params{iParam};

    x = zig_s2p.(p);
    y = zig_mask.(p);

    valid = isfinite(x) & isfinite(y);

    nexttile;
    hold on;

    scatter(x(valid), y(valid), 12, 'filled', ...
        'MarkerFaceAlpha', 0.5);

    mn = min([x(valid); y(valid)]);
    mx = max([x(valid); y(valid)]);

    plot([mn mx], [mn mx], 'k--');

    axis square;
    xlim([mn mx]);
    ylim([mn mx]);

    xlabel(['suite2p ' labels{iParam}]);
    ylabel(['maskNMF ' labels{iParam}]);

    if sum(valid) >= 3
        [rho,pval] = corr(x(valid), y(valid), ...
            'Type','Spearman', ...
            'Rows','complete');

        title(sprintf('\\rho = %.2f, p = %.1g', rho, pval));
    end

    box off;
end

%% plot repeat 1 and repeat 2 of some reliable neurons
% Fixed axes layout shared by both figures
figPos = [100 100 2400 1100];   % pixels: [left bottom width height]
left   = 0.05;
right  = 0.04;
bottom = 0.06;
top    = 0.1;
gap    = 0.006;
axW = 1 - left - right;
axH = (1 - top - bottom - (k-1)*gap) / k;
axPos = zeros(k,4);
for i = 1:k
    axBottom = 1 - top - i*axH - (i-1)*gap;
    axPos(i,:) = [left axBottom axW axH];
end

figure('Position', figPos);
ax = gobjects(k,1);

for i = 1:k
    ax(i) = axes('Position', axPos(i,:));
    hold on;
    y1 = results(iSess).masknmf.traces.repeat1(:,iROI_maskNMF(i));
    y2 = 0.2*max(y1)+results(iSess).masknmf.traces.repeat2(:,iROI_maskNMF(i));
    plot(results(iSess).time,y1,'Color',[.1 .1 .1]);
    plot(results(iSess).time,y2,'Color',[.4 .4 .4]);
    if i<k
        set(gca,'XTick',[])
        ax(i).XAxis.Color = 'none';
    else
        xlabel('time (s)');
    end
    ylims = [0,1.01*max([y1;y2])];
    ylim(ylims);
    set(gca,'YTick',floor(ylims));
    text(300,0.5*ylims(2),sprintf('r=%.2f',bestCorr(iROI_maskNMF(i))),'HorizontalAlignment','left');
    set(ax(i), 'Position', axPos(i,:));  % enforce exact same box
end
linkaxes(ax,'x');
sgtitle(sprintf('maskNMF ROIs (top %d r)',k))

%suite2p
figure('Position', figPos);
ax = gobjects(k,1);

for i = 1:k
    ax(i) = axes('Position', axPos(i,:));
    hold on;
    y1 = results(iSess).s2p.traces.repeat1(:,iROI_s2p(i));
    y2 = 0.2*max(y1)+results(iSess).s2p.traces.repeat2(:,iROI_s2p(i));
    plot(results(iSess).time,y1,'Color',[.1 .1 .1]);
    plot(results(iSess).time,y2,'Color',[.4 .4 .4]);
    if i<k
        set(gca,'XTick',[])
        ax(i).XAxis.Color = 'none';
    else
        xlabel('time (s)');
    end
    ylims = [0,1.01*max([y1;y2])];
    ylim(ylims);
    set(gca,'YTick',floor(ylims));
    text(300,0.5*ylims(2),sprintf('r=%.2f',bestCorr(iROI_maskNMF(i))),'HorizontalAlignment','left');
    set(ax(i), 'Position', axPos(i,:));  % enforce exact same box
end
linkaxes(ax,'x');
sgtitle('suite2p ROIs (putative matches)')

%% Scatter plot comparing reliability of matched ROIs


x = results(iSess).masknmf.r(iMatches);
y = results(iSess).s2p.r(bestS2pIdx_matches);

cMin = 0.5;
cMax = 1.0;

figure;

% Left: scatter plot
ax1 = subplot(1,2,1);

scatter(x, y, 10, matchCorr_matches, 'filled', ...
    'MarkerFaceAlpha', 0.7);

hold on;
plot([-.2 1], [-.2 1], 'k');

axis square;
xlim([-.2 1]);
ylim([-.2 1]);

xlabel('maskNMF cross-repeat reliability');
ylabel('suite2p cross-repeat reliability');
%title('Matched ROIs');

nCmap = 256;
baseMap = parula(nCmap);
fadeFrac = linspace(0,1,nCmap)';
startColor = 0.92 * [1 1 1];
cmap = (1-fadeFrac).*startColor + fadeFrac.*baseMap;

colormap(ax1, cmap);
clim([cMin cMax]);

cb = colorbar;
cb.Label.String = 'Trace-match correlation';

% Right: paired bar plot by trace-match correlation bin
ax2 = subplot(1,2,2);
hold on;

% Bin definitions
binEdges = 0.6:0.05:1.0;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
binLabels = strings(1, numel(binEdges)-1);

pVals = nan(numel(binCenters),1);

for iBin = 1:numel(binCenters)

    lo = binEdges(iBin);
    hi = binEdges(iBin+1);

    if iBin < numel(binCenters)
        inBin = matchCorr_matches >= lo & matchCorr_matches < hi;
    else
        inBin = matchCorr_matches >= lo & matchCorr_matches <= hi;
    end

    binLabels(iBin) = sprintf('%.2f-%.2f', lo, hi);

    xBin = x(inBin);
    yBin = y(inBin);

    if isempty(xBin)
        continue
    end

    avgMask = median(xBin,'omitnan');
    avgS2p  = median(yBin,'omitnan');

    % Color according to trace-match correlation bin
    thisColor = interp1( ...
        linspace(cMin,cMax,size(cmap,1)), ...
        cmap, ...
        binCenters(iBin), ...
        'linear', ...
        'extrap');

    % Individual paired ROIs
    xPair = [iBin-0.25 iBin+0.25];

    for j = 1:numel(xBin)
        plot(xPair,...
            [xBin(j) yBin(j)],...
            '-',...
            'Color',[0 0 0 0.05],...
            'LineWidth',0.25);
    end

    % Thick mean line
    plot(xPair,...
        [avgMask avgS2p],...
        '-',...
        'Color',thisColor,...
        'LineWidth',3);

    % Wilcoxon signed-rank test
    valid = ~isnan(xBin) & ~isnan(yBin);

    if sum(valid) >= 2
        pVals(iBin) = signrank(xBin(valid),yBin(valid));
    end
end

%Significance stars
yl = ylim;
yRange = range(yl);

for iBin = 1:numel(pVals)

    if isnan(pVals(iBin))
        continue
    elseif pVals(iBin) < 0.001
        stars = '***';
    elseif pVals(iBin) < 0.01
        stars = '**';
    elseif pVals(iBin) < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end

    if ~isempty(stars)
        text(iBin,...
            yl(2) + 0.03*yRange,...
            stars,...
            'HorizontalAlignment','center',...
            'FontSize',10,...
            'FontWeight','bold');
    end
end

ylim([yl(1), yl(2)+0.1*yRange]);

yText = yl(2) - 0.2*range(yl);
text(1-0.2, yText, 'maskNMF', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90, ...
    'FontSize', 8);
text(1+0.2, yText, 'suite2p', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90, ...
    'FontSize', 8);

set(gca,...
    'XTick',1:numel(binCenters),...
    'XTickLabel', binLabels);
xtickangle(45)

xlabel('Trace-match correlation');
ylabel('Cross-repeat reliability');
box off;
axis square;

sgtitle('maskNMF v. suite2P reliability, binned by match correlation');


%% shuffle control
for iSess = 1:numel(results)

    [results(iSess).masknmf.rNearestShuffle, ...
     results(iSess).masknmf.nearestIdx, ...
     results(iSess).masknmf.nearestDist] = nearestNeighborRepeatReliability( ...
        results(iSess).masknmf.traces.repeat1, ...
        results(iSess).masknmf.traces.repeat2, ...
        results(iSess).masknmf.fovs, ...
        results(iSess).masknmf.pos);

    [results(iSess).s2p.rNearestShuffle, ...
     results(iSess).s2p.nearestIdx, ...
     results(iSess).s2p.nearestDist] = nearestNeighborRepeatReliability( ...
        results(iSess).s2p.traces.repeat1, ...
        results(iSess).s2p.traces.repeat2, ...
        results(iSess).s2p.fovs, ...
        results(iSess).s2p.pos);

end

%% plot for example sess
figure;
cols = colororder;

subplot(1,2,1);
hold on
histogram(results(iSess).masknmf.r,-.1:.025:0.9,...
    'DisplayStyle','stairs','EdgeColor',cols(2,:),...
    'LineWidth',2,'EdgeAlpha',0.5);
histogram(results(iSess).masknmf.rNearestShuffle,-.1:.025:0.9,...
    'DisplayStyle','stairs','EdgeColor',[.5 .5 .5],...
    'LineWidth',2,'EdgeAlpha',0.5);
ylim([0,2000]);
legend({'true same ROI','nearest-neighbor shuffle'});
xlabel('Cross-repeat reliability');
ylabel('Count');
title('maskNMF');

subplot(1,2,2);
hold on
histogram(results(iSess).s2p.r,-.1:.025:0.9,...
    'DisplayStyle','stairs','EdgeColor',cols(1,:),...
    'LineWidth',2,'EdgeAlpha',0.5);
histogram(results(iSess).s2p.rNearestShuffle,-.1:.025:0.9,...
    'DisplayStyle','stairs','EdgeColor',[.5 .5 .5],...
    'LineWidth',2,'EdgeAlpha',0.5);
ylim([0,2000]);
legend({'true same ROI','nearest-neighbor shuffle'});
xlabel('Cross-repeat reliability');
ylabel('Count');
title('suite2p (matching total nr of ROIs)');

%% plot some traces from shuffled pairs

method = 'masknmf';   % or 's2p'
iSess = 1;
[~,iROI] = max(results(iSess).masknmf.rNearestShuffle);

iNN = results(iSess).(method).nearestIdx(iROI);
iNN = round(iNN);

t = results(iSess).time;

roiRep1 = results(iSess).(method).traces.repeat1(:,iROI);
roiRep2 = results(iSess).(method).traces.repeat2(:,iROI);
nnRep1 = results(iSess).(method).traces.repeat1(:,iNN);
nnRep2  = results(iSess).(method).traces.repeat2(:,iNN);

rs = corrcoef(roiRep1,nnRep1);
r_rep1 = rs(2,1);

figure;

ax(1) = subplot(4,1,1);
plot(t, roiRep1, 'k', 'LineWidth',1.5);
box off

title(ax(1), sprintf('ROI %d repeat 1', iROI));

ax(2) = subplot(4,1,2);
plot(t, roiRep2, 'Color', [.5 .5 .5], 'LineWidth',1.5);
ylabel('deconv. act.');
box off

ax(3) = subplot(4,1,3);
plot(t, nnRep1, 'k', 'LineWidth',1.5);
box off

ax(4) = subplot(4,1,4);
plot(t, nnRep2, 'Color', [.5 .5 .5], 'LineWidth',1.5);
xlabel('time (s)');
box off

linkaxes(ax,'xy');

% Only show x-axis on bottom subplot
for i = 1:3
    ax(i).XAxis.Color = 'none';
    ax(i).XTick = [];
end

ax(4).XAxis.Color = [0 0 0];

% Optional: make y-axis limits identical
yl = [-0.3, ...
      1.02*max([roiRep1; roiRep2; nnRep1; nnRep2])];

ylim(ax(1),yl);
ylim(ax(2),yl);
ylim(ax(3),yl);
ylim(ax(4),yl);

title(ax(2), sprintf('ROI %d repeat 2 (r = %.2f)', ...
    iROI, results(iSess).(method).r(iROI)));

title(ax(3), sprintf('NN ROI %d (d=%.1f px) repeat 1 (r = %.2f)', ...
    iNN, results(iSess).(method).nearestDist(iROI), ...
    r_rep1));

title(ax(4), sprintf('NN ROI %d (d=%.1f px) repeat 2 (r = %.2f)', ...
    iNN, results(iSess).(method).nearestDist(iROI), ...
    results(iSess).(method).rNearestShuffle(iROI)));
