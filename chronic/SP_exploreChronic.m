%root = 'Y:\Subjects\SP061';
root_chronic = 'G:\Shared drives\WG-Mesoscope\Chronic\SP061';
%root_chronic = 'G:\Shared drives\WG-Mesoscope\Sample data\SP058';
%iGoodPaths = [4,5,6,7];
iGoodPaths = [1,2,3,4,5,6,7];
%iGoodPaths = 1:10;

% paths = IBL_listSessionPaths('root',root,'protocol',{'cuedBiasedChoiceWorld'},'mpci',true);
% paths_new = {};
% for i = 1:length(paths)-4 %hack to exclude last one
%     if ~isempty(dir(fullfile(paths{i}, '**', '_sp_video.times.npy')))
%         paths_new{end+1} = paths{i};
%     end
% end

%% find all the ROICaT results files and load the relevant fields
results_fn = '.ROICaT.tracking.results.mat';
roicat_results = dir(fullfile(root_chronic,['*' results_fn]));

fov_names = {};
labels_bySession_all = {};
ROIs_bySession_all = {};
paths_stat_all = {};
for i = 1:length(roicat_results)
    str = roicat_results(i).name;
    idx = strfind(str, results_fn);
    fov_names{end+1} = str(1:idx-1); %make sure we keep track of FOV names (THIS ASSUMES THE FOV NAMES ARE THE SAME ACROSS DAYS!!)
    load(fullfile(roicat_results(i).folder,roicat_results(i).name));
    labels_bySession_all{end+1} = clusters.labels_bySession(iGoodPaths);
    ROIs_bySession_all{end+1} = ROIs.ROIs_aligned(iGoodPaths);
    paths_stat_all = input_data.paths_stat(iGoodPaths);
end

%% now load the traces, neuralFrameTimes and videoFrameTimes of each FOV
cell_template = cell(length(roicat_results),length(labels_bySession_all{1}));
[F_all, neuralFrameTimes_all, videoFrameTimes_all, frameQC_all, badframes_all] = deal(cell_template);
for i = 1:length(roicat_results)
    for j = 1:length(labels_bySession_all{i})
        FOV_dir = dir(fullfile(paths_stat_all{i,j},'..','..','*times.npy'));
        video_dir = dir(fullfile(FOV_dir.folder,'..','**','_sp_video.times.npy'));
        %F = readNPY(fullfile(FOV_dir.folder,'mpci.ROIActivityF.npy'));
        %Fneu = readNPY(fullfile(FOV_dir.folder,'mpci.ROINeuropilActivityF.npy'));
        %F_all{i,j} = F - 0.7*Fneu;
        F_all{i,j} = readNPY(fullfile(FOV_dir.folder,'mpci.ROIActivityDeconvolved.npy'));
        frameQC_all{i,j} = readNPY(fullfile(FOV_dir.folder,'mpci.mpciFrameQC.npy'));
        badframes_all{i,j} = readNPY(fullfile(FOV_dir.folder,'mpci.badFrames.npy'));
        neuralFrameTimes_all{i,j} = readNPY(fullfile(FOV_dir.folder,'mpci.times.npy'));
        videoFrameTimes_all{i,j} = readNPY(fullfile(video_dir(1).folder,video_dir(1).name));
    end
end

%% extract relevant neural traces (at video frames and for full trace)

iFOV = 1; %select one FOV (TODO batch run across all FOVs from a given session)
tau_video = 15; %smoothing window for video (in video frames, relative to lowest frame-rate, usually 30Hz)
tau_full = 1; %smoothing window for full trace (in neural frames)
ds_factor = 10; %downsampling factor (for correlations)
H = 512; %height in pixels
W = 512; %width in pixels

%get all putative cluster idxs that are present across all days
nsessions = length(labels_bySession_all{iFOV});
labels_full = horzcat(labels_bySession_all{iFOV}{:});
[labels_unique,~,idx] = unique(labels_full);
labels_counts = accumarray(idx, 1);
goodLabels_unique = unique(labels_unique(labels_counts == nsessions));
nGoodLabels = length(goodLabels_unique);

%get cluster quality metrics (ROICaT)
cluster_silhouette = [quality_metrics.cluster_silhouette{:}];
cluster_intra_means = [quality_metrics.cluster_intra_means{:}];
confidence = ((cluster_silhouette + 1) / 2) .* cluster_intra_means;
confidence_goodLabels = confidence(goodLabels_unique+2)';

%get minimum video frame length & time (not all videos are the same frameRate)
[nVideoFrames_min, imin] = min(cellfun(@length,videoFrameTimes_all));
tVideo = videoFrameTimes_all{imin}(:,1) - videoFrameTimes_all{imin}(1,1);

%check if there are badframes anywhere in the video
for iSess = 1:nsessions
    bad_frames_bool = (badframes_all{iFOV,iSess} | frameQC_all{iFOV,iSess}>0);
    bad_frames_times = neuralFrameTimes_all{iFOV,iSess}(bad_frames_bool);
    video_start_end_times = [videoFrameTimes_all{iFOV,iSess}(1,1) - videoFrameTimes_all{iFOV,iSess}(end,2)];
    if any(bad_frames_bool)
        if any(bad_frames_times>video_start_end_times(1) && bad_frames_times<video_start_end_times(2))
            warning('There appear to be bad frames in the video times, be careful!')
        end
    end
end

%resample neural trace at each video frame (in each repeat)
F_video = nan(nGoodLabels,nsessions,nVideoFrames_min,2);
ROI_masks = {};
cluster_centroids = {};
F_video_smooth = nan(nGoodLabels,nsessions,nVideoFrames_min/ds_factor,2);
for iCluster=1:nGoodLabels
    for iSess = 1:nsessions
        iROI = labels_bySession_all{iFOV}{iSess}==goodLabels_unique(iCluster);
        for iRep = 1:2
            nVideoFrames = size(videoFrameTimes_all{iFOV,iSess},1);
            frRatio = nVideoFrames/nVideoFrames_min;
            videoFrameTimes = videoFrameTimes_all{iFOV,iSess}(round(linspace(1,(nVideoFrames-ceil(frRatio/2)),nVideoFrames_min)),iRep);
            F = interp1(neuralFrameTimes_all{iFOV,iSess}, F_all{iFOV,iSess}(:,iROI), videoFrameTimes, 'linear');
            F_video(iCluster,iSess,:,iRep) = F;
            F_video_smooth(iCluster,iSess,:,iRep) = resample(smooth(double(F),tau_video),1,ds_factor);
        end
    end
end

%get smoothed & downsampled activity of the reference population (excluding the bad frames)
F_full_smooth = {};
for iSess = 1:nsessions
    frames_beforeVideo = neuralFrameTimes_all{iFOV,iSess} < (videoFrameTimes_all{iFOV,iSess}(1,1) - 10); %assume passiveVideo block started 10 seconds before first video frame
    frames_toinclude = ((~badframes_all{iFOV,iSess} | frameQC_all{iFOV,iSess}==0) & frames_beforeVideo); %exclude bad frames & passiveVideo block
    F = F_all{iFOV,iSess}(frames_toinclude,:); 
    F_full_smooth{iSess} = nan(nGoodLabels,ceil(sum(frames_toinclude)/ds_factor));
    for iCluster=1:nGoodLabels
        iROI = labels_bySession_all{iFOV}{iSess}==goodLabels_unique(iCluster);
        F_full_smooth{iSess}(iCluster,:) = resample(smooth(double(F(:,iROI)),tau_full),1,ds_factor);
    end
end

%% compute nearest neighbour for each cluster

%sub-select good clusters from ROIs_bySession_all and order by cluster ID
clusters_bySession_all = {{}};
for iSess = 1:nsessions
    for iCluster=1:nGoodLabels
        iROI = labels_bySession_all{iFOV}{iSess}==goodLabels_unique(iCluster);
        clusters_bySession_all{iFOV}{iSess}(iCluster,:) = ROIs_bySession_all{iFOV}{iSess}(iROI,:);
    end
end

%find nearest neighbour of each cluster
iNearestCluster = {};
for iSess = 1:nsessions
        sparse_data = clusters_bySession_all{iFOV}{iSess};
        iNearestCluster{iSess} = find_nearest_neighbours(sparse_data);
end


%% for each cluster, compute correlation between video repeats (either with itself or its nearest neighbour)

corrs_same = nan(nGoodLabels,nsessions,nsessions);
corrs_diff = nan(nGoodLabels,nsessions,nsessions);
for i = 1:nGoodLabels
    for iSess = 1:nsessions
        for jSess = 1:nsessions
            corrs_same(i,iSess,jSess) = corr(squeeze(F_video_smooth(i,iSess,:,1)),squeeze(F_video_smooth(i,jSess,:,2)));
            corrs_diff(i,iSess,jSess) = corr(squeeze(F_video_smooth(i,iSess,:,1)),squeeze(F_video_smooth(iNearestCluster{iSess}(i),jSess,:,2)));
        end
    end
end

corrs_within = nan(nGoodLabels,nsessions);
corrs_within_mean = nan(nGoodLabels,1);
corrs_between_mean = nan(nGoodLabels,1);
corrs_within_each = nan(nGoodLabels,nsessions);
corrs_between_each = nan(nGoodLabels,nsessions);
corrs_within_diff = nan(nGoodLabels,nsessions);
for iSess = 1:nsessions
    corrs_within(:,iSess)=corrs_same(:,iSess,iSess);
end
I = logical(eye(nsessions));
for iCluster = 1:nGoodLabels
    corrs = squeeze(corrs_same(iCluster,:,:));
    corrs_d = squeeze(corrs_diff(iCluster,:,:));
    corrs_within_mean(iCluster) = mean(corrs(I),'all');
    corrs_between_mean(iCluster) = mean(corrs(~I),'all');
    corrs_within_each(iCluster,:) = corrs(I);
    corrs_within_diff(iCluster,:) = corrs_d(I);
    for iSess = 1:nsessions
        idxs = true(1,nsessions);
        idxs(iSess) = false;
        corrs_between_each(iCluster,iSess) = mean([corrs(iSess,idxs),corrs(idxs,iSess)']);
    end
end

%% for each cluster, compute instantaneous cross-correlation against the reference population (across full recording)
iRefPop = find(confidence_goodLabels>0.75);
nRefPop = length(iRefPop);
crosscorrs = nan(nsessions,nGoodLabels,nRefPop,2);
for iSess = 1:nsessions
    nF = size(F_full_smooth{iSess},2);
    tic
    for i=1:nGoodLabels
        for j=1:nRefPop
            if i~=iRefPop(j)
                crosscorrs(iSess,i,j,1) = corr(squeeze(F_full_smooth{iSess}(i,1:floor(nF/2)))',squeeze(F_full_smooth{iSess}(iRefPop(j),1:floor(nF/2)))');
                crosscorrs(iSess,i,j,2) = corr(squeeze(F_full_smooth{iSess}(i,ceil(nF/2):end))',squeeze(F_full_smooth{iSess}(iRefPop(j),ceil(nF/2):end))');
            end
        end
    end
    toc
end

refcorr_same = nan(nGoodLabels,nsessions,nsessions);
refcorr_diff = nan(nGoodLabels,nsessions,nsessions);
for i = 1:nGoodLabels
    for iSess = 1:nsessions
        for jSess = 1:nsessions
            refcorr_same(i,iSess,jSess) = corr(squeeze(crosscorrs(iSess,i,:,1)),squeeze(crosscorrs(jSess,i,:,2)),'rows','complete');
            refcorr_diff(i,iSess,jSess) = corr(squeeze(crosscorrs(iSess,i,:,1)),squeeze(crosscorrs(jSess,iNearestCluster{iSess}(i),:,2)),'rows','complete');
        end
    end
end

refcorr_within = nan(nGoodLabels,nsessions);
refcorr_within_mean = nan(nGoodLabels,1);
refcorr_between_mean = nan(nGoodLabels,1);
refcorr_within_each = nan(nGoodLabels,nsessions);
refcorr_between_each = nan(nGoodLabels,nsessions);
refcorr_within_diff = nan(nGoodLabels,nsessions);
for iSess = 1:nsessions
    refcorr_within(:,iSess)=refcorr_same(:,iSess,iSess);
end
I = logical(eye(nsessions));
for iCluster = 1:nGoodLabels
    refcorr = squeeze(refcorr_same(iCluster,:,:));
    refcorr_d = squeeze(refcorr_diff(iCluster,:,:));
    refcorr_within_mean(iCluster) = mean(refcorr(I),'all');
    refcorr_between_mean(iCluster) = mean(refcorr(~I),'all');
    refcorr_within_each(iCluster,:) = refcorr(I);
    refcorr_within_diff(iCluster,:) = refcorr_d(I);
    for iSess = 1:nsessions
        idxs = true(1,nsessions);
        idxs(iSess) = false;
        refcorr_between_each(iCluster,iSess) = mean([refcorr(iSess,idxs),refcorr(idxs,iSess)']);
    end
end

%% for each cluster, correlate its activity on repeat 1 with activity from all neurons on repeat 2 (within and across days)
%...and compute the percentiles/ranks of within-neuron correlation relative to cross-neuron correlations

corr_prctiles = nan(nGoodLabels,nsessions,nsessions);
corr_ranks = nan(nGoodLabels,nsessions,nsessions);
corr_full_r  = nan(nGoodLabels,nGoodLabels,nsessions,nsessions);
corr_full_p  = nan(nGoodLabels,nGoodLabels,nsessions,nsessions);
for i = 1:nGoodLabels
    if mod(i,10)==0
        i
    end
    for iSess = 1:nsessions
        for jSess = 1:nsessions
            cnt = 0;
            corr_x = nan(1,nGoodLabels);
            for j = 1:nGoodLabels
                [x,p] = corr(squeeze(F_video_smooth(i,iSess,:,1)),squeeze(F_video_smooth(j,jSess,:,2)),'type','Spearman','tail','right');
                if i==j
                    %corrs_same(i,iSess,jSess) = corr(smooth(squeeze(F_video(i,iSess,:,1)),tau),smooth(squeeze(F_video(j,jSess,:,2)),tau));
                    corr_0 = x;
                else
                    cnt = cnt+1;
                    %corrs_diffs(i,iSess,jSess,cnt) = corr(smooth(squeeze(F_video(i,iSess,:,1)),tau),smooth(squeeze(F_video(j,jSess,:,2)),tau));
                    corr_x(cnt) = x;
                end
                %if iSess==1
                corr_full_r(i,j,iSess,jSess) = x;
                corr_full_p(i,j,iSess,jSess) = p;
                %end
            end
            corr_prctiles(i,iSess,jSess) = comp_percentile(corr_x,corr_0);
            rank = tiedrank([corr_0,corr_x]);
            corr_ranks(i,iSess,jSess) = nGoodLabels-rank(1)+1;
        end
    end
end

%% select visually reliable neurons (within day)
corr_ps_all_within = nan(nGoodLabels,nsessions);
corr_rs_all_within = nan(nGoodLabels,nsessions);
corr_ranks_all_within = nan(nGoodLabels,nsessions);
corr_prctiles_all_within = nan(nGoodLabels,nsessions);
I = logical(eye(nsessions)); %identity matrix
for i = 1:nGoodLabels
    c_r = squeeze(corr_full_r(i,i,:,:));
    c_rank = squeeze(corr_ranks(i,:,:));
    c_p = squeeze(corr_full_p(i,i,:,:));
    c_ptile = squeeze(corr_prctiles(i,:,:));
    corr_prctiles_all_within(i,:) = c_ptile(I);
    corr_ranks_all_within(i,:) = c_rank(I);
    corr_ps_all_within(i,:) = c_p(I);
    corr_rs_all_within(i,:) = c_r(I);
end

%select day
D = 1;

%select criteria

%iGood = corr_rs_all_within(:,D)>0.1 & corr_ranks_all_within(:,D)<=1;
%iGood = corr_rs_all_within(:,D)>0 & corr_prctiles_all_within(:,D)>95;

% iGood = corr_rs_all_within(:,D)>0.1; %neurons that had large correlation coefficient with themselves on day D (r>0.2) 
% iGood = corr_ps_all_within(:,D)<0.005; %neurons that were significantly correlated with themselves on day D (p<0.005)
% iGood = any(corr_prctiles_all_within>90); %neurons that had large correlation with themselves (compared to other neurons) on any day
% iGood = corr_prctiles_all_within(:,D)>90; %neurons that had large correlation with themselves (compared to other neurons) on day D
% iGood = any(corr_ranks_all_within<=1,2); %neurons that were most correlated to themselves on any day
iGood = corr_ranks_all_within(:,D)<=1; %neurons that were most correlated to themselves on day D

% OPTIONAL: additionally select based on a quality metric of the cluster
hiConfidence = confidence_goodLabels>0.5;
%iGood = iGood & hiConfidence;
sum(iGood)

%% for an example neuron, plot the traces against each other

%find 'typical example neuron'
iThis = D; %choose day
%[m,iExampleROI] = max(squeeze(corr_rs_all_within(:,iThis))); %neuron with largest corr coef on this day
[m,i] = max(squeeze(corr_rs_all_within(iGood,iThis))); goodROIs = find(iGood); iExampleROI = goodROIs(i); %'reliable' neuron with largest corr coef

%[m,iTop] = max(corrs_within_mean);

%iExampleROI = 276;
iNNeighbour = mode(cellfun(@(x) x(iExampleROI), iNearestCluster));


%plot this neuron's traces (within and across days)

%maxval = max([F_this_1;F_this_2;F_diff_1;F_diff_2]);

figure;
cols = colororder;
t = tiledlayout(nsessions,9);
ax = [];
axx = [];
for i = 1:nsessions
    ax(i)=nexttile(t,[1,9]);
    hold on;
    plot(tVideo,squeeze(F_video(iExampleROI,i,:,1)),'color',cols(1,:));
    plot(tVideo,squeeze(F_video(iExampleROI,i,:,2)),'color',cols(2,:));
    ylims = get(gca,'ylim'); maxval = ylims(2);
    text(10,0.9*maxval,sprintf('r=%.3f',corrs_same(iExampleROI,i,i)));
    %axx(i) = nexttile(t,[1,1]);
    %imagesc(squeeze(ROI_masks{iTop}(i,:,:)));
    %axis off
end
linkaxes(ax,'x');
%linkaxes(ax,'xy');

%nexttile;
%hold on;
%F_this_1 = squeeze(F_video(iExampleROI,iThis,:,1));
%F_diff_2 = squeeze(F_video(iExampleROI,iDiff,:,2));
%plot(tVideo,F_this_1,'color',cols(1,:));
%plot(tVideo,F_diff_2,'color',cols(4,:));
%text(10,0.9*maxval,sprintf('r=%.3f',corrs_tolast(iTop,iDiff)));

%plot correlation coefficient matrix of this neuron
figure('Position',[100,100,300,400]);
tiledlayout(2,1);
corrMat = squeeze(corrs_same(iExampleROI,:,:));
corrMat_diff = squeeze(corrs_diff(iExampleROI,:,:));
nexttile;
imagesc(corrMat);
caxis([0,max(corrMat(:))]);
colorbar;
xlabel('day (repeat 1)');
ylabel('day (repeat 2)');
title(sprintf('cluster %d v cluster %d',iExampleROI,iExampleROI));
nexttile;
imagesc(corrMat_diff);
caxis([0,max(corrMat(:))]);
colorbar;
xlabel('day (repeat 1)');
ylabel('day (repeat 2)');
%title(sprintf('cluster %d v cluster %d',iTop,iNNeighbour));
title(sprintf('cluster %d v nearest neighbour (%d)',iExampleROI,iNNeighbour));

%% across all clusters, plot within v. across day correlation coefficients

%iGood = any(corrs_within>0.4,2);

figure('Units','normalized','Position',[0.26 0.36 0.1 0.5]);
tiledlayout(nsessions,1);
for iSess = 1:nsessions
    nexttile;
    hold on
    scatter(corrs_within_each(iGood,iSess),corrs_between_each(iGood,iSess));
    scatter(corrs_within_each(iExampleROI,iSess),corrs_between_each(iExampleROI,iSess),'filled');
    plot([-.1,1],[-.1,1],'-k')
    xlim([-.1,1]);
    ylim([-.1,1]);
    axis square
end

figure('Units','normalized','Position',[0.36 0.36 0.1 0.26]);
tiledlayout(3,1);

nexttile([1,1]);
hold on;
c1 = corrs_within_each(iGood,:); plot_histogram_line(c1(:),[-1:.025:1]);
c2 = corrs_between_each(iGood,:); plot_histogram_line(c2(:),[-1:.025:1]);
c3 = corrs_within_diff(iGood,:); plot_histogram_line(c3(:),[-1:.025:1]);
ylims = get(gca,'ylim');
xlims = get(gca,'xlim');
text(xlims(1)+0.05*(diff(xlims)),ylims(1)+0.85*(diff(ylims)),'same, within Ds','Color',cols(1,:))
text(xlims(1)+0.05*(diff(xlims)),ylims(1)+0.65*(diff(ylims)),'matched, across Ds','Color',cols(2,:))
text(xlims(1)+0.05*(diff(xlims)),ylims(1)+0.75*(diff(ylims)),'different, within Ds','Color',cols(3,:))
xlabel('correlation')
ylabel('proportion')
%histogram(corrs_within_each(:),[-.1:.025:0.9],'Normalization','probability','DisplayStyle','stairs');
%histogram(corrs_between_each(:),[-.1:.025:0.9],'Normalization','probability','DisplayStyle','stairs');
%histogram(corrs_within_diff(:),[-.1:.025:0.9],'Normalization','probability','DisplayStyle','stairs');
%
%hold on
%scatter(corrs_within_mean(iGood),corrs_between_mean(iGood));
%plot([-.1,.5],[-.1,.5],'-k')
%scatter(corrs_within
%axis square

%plot ROC curves and get AUCs of 
% (1) same-day-same-ROIs against same-day-different-ROIs 
% (2) different-day-matched-ROIs against same-day-different-ROIs
[X_within, Y_within, ~, AUC_within] = perfcurve([zeros(length(c3(:)), 1); ones(length(c1(:)), 1)], [c3(:); c1(:)], 1);
[X_between, Y_between, ~, AUC_between] = perfcurve([zeros(length(c3(:)), 1); ones(length(c2(:)), 1)], [c3(:); c2(:)], 1);
nexttile([2,1]);
plot(X_within, Y_within, 'LineWidth', 2); hold on;
plot(X_between, Y_between, 'LineWidth', 2);
plot([0, 1], [0, 1], 'k--');
text(0.7,0.35,sprintf('AUC=%.3f',AUC_within),'Color',cols(1,:))
text(0.7,0.3,sprintf('AUC=%.3f',AUC_between),'Color',cols(2,:))
axis square
xlabel('False positives')
ylabel('Hits')

%disp(['AUC_within: ', num2str(AUC_within)]);
%disp(['AUC_between: ', num2str(AUC_between)]);

%% across all clusters, plot results of cross correlations with reference population 

figure('Units','normalized','Position',[0.3 0.36 0.5 0.2]);
tiledlayout(1,7);
cmap = brewermap([],'*RdYlBu');
ax = [];
for iSess = 1:nsessions
    ax(iSess) = nexttile;
    imagesc(squeeze(crosscorrs(iSess,confidence_goodLabels>0.75,:,1)));
    colormap(cmap)
    caxis([-.3,.3])
    %set(gca,'XTick',[1,nGoodLabels],'YTick',[1,nGoodLabels])
    set(gca,'XTick',[1,50],'YTick',[1,50])
    %caxis([0,0.6]);
    axis square
    if iSess==1
        xlabel('cluster (first half)')
        ylabel('cluster (second half)')
    end
    title(sprintf('Day %d',iSess))
    if iSess == nsessions
        cb = colorbar;
        cb.Label.String = 'r';
    end
end
linkaxes(ax,'xy');

figure('Units','normalized','Position',[0.36 0.36 0.1 0.26]);
tiledlayout(3,1);

hiConfidence = confidence_goodLabels>0.5;

nexttile([1,1]);
hold on;
c1 = refcorr_within_each(hiConfidence,:); plot_histogram_line(c1(:),[-1:.025:1]);
c2 = refcorr_between_each(hiConfidence,:); plot_histogram_line(c2(:),[-1:.025:1]);
c3 = refcorr_within_diff(hiConfidence,:); plot_histogram_line(c3(:),[-1:.025:1]);
ylims = get(gca,'ylim');
xlims = get(gca,'xlim');
text(xlims(1)+0.05*(diff(xlims)),ylims(1)+0.85*(diff(ylims)),'same, within Ds','Color',cols(1,:))
text(xlims(1)+0.05*(diff(xlims)),ylims(1)+0.65*(diff(ylims)),'matched, across Ds','Color',cols(2,:))
text(xlims(1)+0.05*(diff(xlims)),ylims(1)+0.75*(diff(ylims)),'different, within Ds','Color',cols(3,:))
xlabel('correlation')
ylabel('proportion')

%plot ROC curves and get AUCs of 
% (1) same-day-same-ROIs against same-day-different-ROIs 
% (2) different-day-matched-ROIs against same-day-different-ROIs
[X_within, Y_within, ~, AUC_within] = perfcurve([zeros(length(c3(:)), 1); ones(length(c1(:)), 1)], [c3(:); c1(:)], 1);
[X_between, Y_between, ~, AUC_between] = perfcurve([zeros(length(c3(:)), 1); ones(length(c2(:)), 1)], [c3(:); c2(:)], 1);
nexttile([2,1]);
plot(X_within, Y_within, 'LineWidth', 2); hold on;
plot(X_between, Y_between, 'LineWidth', 2);
plot([0, 1], [0, 1], 'k--');
text(0.7,0.35,sprintf('AUC=%.3f',AUC_within),'Color',cols(1,:))
text(0.7,0.3,sprintf('AUC=%.3f',AUC_between),'Color',cols(2,:))
axis square
xlabel('False positives')
ylabel('Hits')

%disp(['AUC_within: ', num2str(AUC_within)]);
%disp(['AUC_between: ', num2str(AUC_between)]);


%% plot results of example neuron

%plot prctiles
figure;
corr_prctiles_iTop = squeeze(corr_prctiles(iExampleROI,:,:));
imagesc(corr_prctiles_iTop); 
caxis([95,100]); 
colormap('gray'); colorbar;

%plot ranks
figure;
range = [1,10];
cmap = flipud(hot(diff(range)));
cmap = [cmap; 0 0 0];
corr_ranks_iTop = squeeze(corr_ranks(iExampleROI,:,:));
imagesc(corr_ranks_iTop); 
caxis([range(1)-0.5, range(end)+0.5]); 
colormap(cmap); cb = colorbar('Direction','reverse');
cb.Label.String = 'Rank';
xlabel('day (repeat 1)')
ylabel('day (repeat 2)')
xticks(1:nsessions)
yticks(1:nsessions)
axis square

%% plot results across neurons

%OPTION 1: full matrix of nsessions x nsessions
% figure('Units','normalized','Position',[0.6 0.36 0.1 0.5]);
% tiledlayout(4,4);
% cmap = brewermap([],'*RdYlBu');
% ax = [];
% for iSess = 1:nsessions
%     for jSess = 1:nsessions
%         ax(iSess,jSess) = nexttile;
%         imagesc(corr_full_r(iGood,iGood,iSess,jSess));
%         colormap(cmap)
%         caxis([-.7,.7])
%         set(gca,'XTick',[],'YTick',[])
%         %caxis([0,0.6]);
%         %colorbar;
%         axis square
%     end
% end
% linkaxes(ax,'xy');

%OPTION 2: session D v. all other sessions
D=1;
%iGood = corr_ranks_all_within(:,D)<=1; %neurons that were most correlated to themselves on day D
figure('Units','normalized','Position',[0.3 0.36 0.3 0.25]);
tiledlayout(2,4,'tilespacing','compact','padding','compact');
cmap = brewermap([],'*RdYlBu');
ax = [];
for iSess = 1:nsessions
    ax(iSess) = nexttile;
    imagesc(corr_full_r(iGood,iGood,iSess,iSess));
    colormap(cmap)
    caxis([-.55,.55])
    set(gca,'XTick',[1,sum(iGood)],'YTick',[1,sum(iGood)])
    %caxis([0,0.6]);
    axis tight square
    if ismember(iSess,[1,5])
        ylabel('cluster (repeat 2)')
    end
    if iSess>4
        xlabel('cluster (repeat 1)')
    end
    title(sprintf('Day %d',iSess))
    if iSess == nsessions
        cb = colorbar;
        cb.Label.String = 'r';
    end
end
linkaxes(ax,'xy');

%%
% corr_prctiles_all_within = nan(1,nGoodLabels);
% for i = 1:nGoodLabels
%     c = squeeze(corr_prctiles(i,:,:));
%     corr_prctiles_max_within(i) = max(c(I));
% end
% iGood = corr_prctiles_max_within>99.9;

corr_prctiles_mean = squeeze(mean(corr_prctiles(iGood,:,:),1));
figure; hold on;
imagesc(corr_prctiles_mean); 
set(gca,'YDir','reverse')
caxis([50,100])
axis tight square
%colormap('gray'); 
colorbar;

%