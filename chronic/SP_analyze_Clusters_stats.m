%this script loads a preprocessed ROI table, with chronically recorded ROIs
%(this should include clusterUIDs). We then do a number of different
%exploratory data analyses using the PETHs and the corresponding test
%statistics


%Define PETH parameters
evnt = 'stimOn'; %event to lock activity to e.g. 'stimOn', 'goCue' or 'feedback'
trialTypeField = 'probabilityLeft'; %define trial type to look at e.g. 'contrastDiff' or 'probabilityLeft' (must match a field in trials table!)
suffix = '_chronicSP067';
stat_to_use = 'priorContrastL';
date_to_load = '250709';
project_name = 'mesoscope_active'; %for loading

pthresh = 0.05;

%define root folders
root = 'Y:\Subjects\'; %server that has the data
root_local = 'C:\Users\Samuel\Documents\2PI\';
saveto = project_name;

%% load the ROI table
ff = fullfile(root_local,project_name,'analysis','PETH');
fn = sprintf('PETH_%s_%s_AllROIs_%s_%s%s.mat',evnt,date_to_load,stat_to_use,trialTypeField,suffix);
load(fullfile(ff,fn));

fn_contrast = sprintf('PETH_%s_%s_AllROIs_slopes_left_contrastDiff%s',evnt,date_to_load,suffix);
s = load(fullfile(ff,fn_contrast),'AllROIs');
AllROIs_contrast = s.AllROIs;

T = linspace(-0.82,3.28,21);

    
%% get cluster peth matrices with associated p-values and test-statistics
[clusterPETHs, clusterPs, clusterVs, clusterPos, clusterLookup, sessionLookup] = get_ClusterMat(AllROIs);
[~, clusterPs_contrast] = get_ClusterMat(AllROIs_contrast);


%% find how reliable the response was to passive video (within day)

% find unique session keys
sessionKeys = arrayfun(@(d) sprintf('%s\\%s\\%s', d.subject, d.date, d.session), AllROIs, 'UniformOutput', false);
uniqueSessions = unique(sessionKeys);

%split paths to get subjects, dates, sessions
splitPaths = split(uniqueSessions,filesep);
if size(splitPaths,2)==1
    subjects = splitPaths(end-2);
    dates = splitPaths(end-1);
    sessions = splitPaths(end);
else
    subjects = splitPaths(:,:,end-2);
    dates = splitPaths(:,:,end-1);
    sessions = splitPaths(:,:,end);
end

% for each cluster in each session, find how reliable the response to the
% passiveVideo was
clusterPassiveR = nan(size(clusterVs));
for iSess = 1:length(sessions)
    
    subject = subjects{iSess};
    date = dates{iSess};
    session = sessions{iSess};
    
    %load the traces
    fprintf('Loading raw 2PI traces...');
    Fall = IBL_loadMesoData_onlyF(subject,date,session);
    
    %load the video frame times
    video_dir = dir(fullfile(root,subject,date,session,'**','_sp_video.times.npy'));
    videoFrameTimes = readNPY(fullfile(video_dir(1).folder,video_dir(1).name));
    
    %get the start and stop frames of the video
    iVideoStartStop = nan(2);
    for i = 1:2
        [~,iVideoStartStop(1,i)] = min(abs(Fall.time-videoFrameTimes(1,i)));
        [~,iVideoStartStop(2,i)] = min(abs(Fall.time-videoFrameTimes(end,i)));
    end
    
    %adjust this in case the two lengths are different
    minLength = min(diff(iVideoStartStop,[],1));
    iVideoStartStop(2,:) = iVideoStartStop(1,:)+minLength;
    
    %convert to logical indices
    iVideoFrameTimes = false(2,size(Fall.spks,2));
    iVideoFrameTimes(1,iVideoStartStop(1,1):iVideoStartStop(2,1)) = true;
    iVideoFrameTimes(2,iVideoStartStop(1,2):iVideoStartStop(2,2)) = true;

    %compute correlation of the traces on two repeats
    for iC = 1:length(clusterLookup)
        iROI = Fall.clusterUID == clusterLookup(iC);
        trace1 = Fall.spks(Fall.clusterUID == clusterLookup(iC),iVideoFrameTimes(1,:));
        trace2 = Fall.spks(Fall.clusterUID == clusterLookup(iC),iVideoFrameTimes(2,:));
        R = corrcoef(trace1,trace2);
        clusterPassiveR(iC,iSess) = R(2,1);
    end

end


%% plot ROI positions

%load allenCCF labels and colors
CCFdir = 'C:\Users\Samuel\Documents\GitHub\allenCCF\';
av = readNPY(fullfile(CCFdir,'annotation_volume_10um_by_index.npy'));
st = loadStructureTree(fullfile(CCFdir,'structure_tree_safe_2017.csv'));
load(fullfile(CCFdir,'Browsing Functions','allen_ccf_colormap_2017.mat'),'cmap');

%load top down atlas
bas = aratopdown.atlas.build_topdown;

%get one example cluster
example_cUID = find(clusterLookup == '1fbdba57-9d29-4f11-abf2-28487c213f24');

%fpos = [50, 100, 1400, 170];
figure('Position',[50, 100 2500 350]);
tiledlayout(1,size(clusterPs,2));

for i=1:size(clusterPs,2)
    
    ax(i) = nexttile;    
    hold on
    
    %first plot atlas boundaries
    cellfun(@(x) cellfun(@(x) plot(1000*x(:,2),1000*x(:,1),'color',[.5 .5 .5]),x,'uni',false), ...
        {bas.dorsal_brain_areas.boundaries_stereotax},'uni', false);
    
    %now plot each ROI in a scatter-plot
    scatter(clusterPos(:,i,1),clusterPos(:,i,2),5,clusterPassiveR(:,i),...
        'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',.3);%,'MarkerEdgeAlpha',.3);
    
    %plot example ROI as a black circle
    scatter(clusterPos(example_cUID,i,1),clusterPos(example_cUID,i,2),10,[0 0 0]);
    
    %colormap(brewermap([],'Reds')); caxis([0,0.5])
    colormap(brewermap([],'*RdYlBu')); caxis([-0.4,0.4]);
    
    axis square
    axis equal
    axis off
    
    %xlim([0,5500]); ylim([-5000,1000]);
    xlim([2700,5400]); ylim([-3700, -1000]);
    daspect([1 1 1])
    
    if i == size(clusterPs,2)
        colorbar;
    end
    
end

linkaxes(ax,'xy');

%% plot visual reliability against beta-weight

figure('Position',fpos); 
for i = 1:7
    subplot(1,7,i);
    hold on;
    scatter(clusterPassiveR(~iSign_pos(:,i) & ~iSign_neg(:,i),i), clusterVs(~iSign_pos(:,i) & ~iSign_neg(:,i),i));
    scatter(clusterPassiveR(iSign_pos(:,i),i), clusterVs(iSign_pos(:,i),i));
    scatter(clusterPassiveR(iSign_neg(:,i),i), clusterVs(iSign_neg(:,i),i));
    %ylim([-1000,1000])
    %xlim([0,1])
end

%% plot histograms of p-values

iSign_pos = clusterPs>1-pthresh;% & clusterVs>0;
iSign_neg = clusterPs<pthresh;% & clusterVs<0;

%iSel = clusterPs_contrast>1-pthresh;
%iSel = clusterPs_contrast<pthresh;

%iSel = clusterPassiveR>0.2;
iSel = false(size(clusterPassiveR));
iSel(mean(clusterPassiveR,2)>0.2,:) = true;
%iSel(median(clusterPassiveR,2)>0.2,:) = true;

cols = colororder;

% plot p-value histograms
fpos = [50, 100, 1400, 170];
figure('Position',fpos); 
for i=1:size(clusterPs,2)
    subplot(1,size(clusterPs,2),i); 
    hold on;
    histogram(clusterPs(~iSign_pos(:,i) & ~iSign_neg(:,i) & iSel(:,i),i),[0:0.025:1],'edgecolor','none'); 
    histogram(clusterPs(iSign_pos(:,i) & iSel(:,i),i),[0:0.025:1],'edgecolor','none'); 
    histogram(clusterPs(iSign_neg(:,i) & iSel(:,i),i),[0:0.025:1],'edgecolor','none'); 
    prop_pos = sum(iSign_pos(:,i) & iSel(:,i))/sum(iSel(:,i));
    prop_neg = sum(iSign_neg(:,i) & iSel(:,i))/sum(iSel(:,i));
    ylim([0,80]);
    ylims = get(gca,'ylim');
    text(0.95,0.8*ylims(2),sprintf('%.1f%%',100*prop_pos),'color',cols(2,:),'HorizontalAlignment','right');
    text(0.05,0.8*ylims(2),sprintf('%.1f%%',100*prop_neg),'color',cols(3,:),'HorizontalAlignment','left');
    hold on; 
    box off; 
    if i>1
        %set(gca,'ytick',[])
        set(gca,'YColor','none'); 
    end
end

%% plot distribution of empirical coefficient


iSign_pos = clusterPs>1-pthresh;% & clusterVs>0;
iSign_neg = clusterPs<pthresh;% & clusterVs<0;

fpos = [50, 100, 1400, 170];
figure('Position',fpos); 
for i=1:size(clusterPs,2)
    subplot(1,size(clusterPs,2),i); 
    hold on;
    histogram(clusterVs(~iSign_pos(:,i) & ~iSign_neg(:,i) & iSel(:,i),i),[-60:4:60],'edgecolor','none'); 
    histogram(clusterVs(iSign_pos(:,i) & iSel(:,i),i),[-60:4:60],'edgecolor','none'); 
    histogram(clusterVs(iSign_neg(:,i) & iSel(:,i),i),[-60:4:60],'edgecolor','none'); 
    %prop_pos = sum(iSign_pos(:,i) & iSel(:,i))/sum(iSel(:,i));
    %prop_neg = sum(iSign_neg(:,i) & iSel(:,i))/sum(iSel(:,i));
    %ylim([0,1080]);
    ylims = get(gca,'ylim');
    %text(0.95,0.8*ylims(2),sprintf('%.1f%%',100*prop_pos),'color',cols(2,:),'HorizontalAlignment','right');
    %text(0.05,0.8*ylims(2),sprintf('%.1f%%',100*prop_neg),'color',cols(3,:),'HorizontalAlignment','left');
    hold on; 
    box off; 
    %set(gca,'YColor','none'); 
end

%% plot ROI positions

%load allenCCF labels and colors
CCFdir = 'C:\Users\Samuel\Documents\GitHub\allenCCF\';
av = readNPY(fullfile(CCFdir,'annotation_volume_10um_by_index.npy'));
st = loadStructureTree(fullfile(CCFdir,'structure_tree_safe_2017.csv'));
load(fullfile(CCFdir,'Browsing Functions','allen_ccf_colormap_2017.mat'),'cmap');

%load top down atlas
bas = aratopdown.atlas.build_topdown;

%fpos = [50, 100, 1400, 170];
figure('Position',[50, 100 2500 350]);
tiledlayout(1,size(clusterPs,2));

%iSel = clusterPs_contrast>-Inf;

for i=1:size(clusterPs,2)
    
    ax(i) = nexttile;    
    hold on
    
    %first plot atlas boundaries
    cellfun(@(x) cellfun(@(x) plot(1000*x(:,2),1000*x(:,1),'color',[.5 .5 .5]),x,'uni',false), ...
        {bas.dorsal_brain_areas.boundaries_stereotax},'uni', false);
    
    %now plot each ROI in a scatter-plot
    %scatter(clusterPos(~iSign_pos(:,i) & ~iSign_neg(:,i),i,1),...
    %    clusterPos(~iSign_pos(:,i) & ~iSign_neg(:,i),i,2),5,[.5 .5 .5],...
    %    'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',.2);%,'MarkerEdgeAlpha',.3);
    scatter(clusterPos(~iSign_pos(:,i) & ~iSign_neg(:,i) & iSel(:,i),i,1),...
        clusterPos(~iSign_pos(:,i) & ~iSign_neg(:,i) & iSel(:,i),i,2),5,[.8,.8,.8],...%cols(1,:),...
        'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',.5);%,'MarkerEdgeAlpha',.3);
    scatter(clusterPos(iSign_pos(:,i) & iSel(:,i),i,1),...
        clusterPos(iSign_pos(:,i) & iSel(:,i),i,2),10,cols(2,:),...
        'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',.5);%,'MarkerEdgeAlpha',.3);
    scatter(clusterPos(iSign_neg(:,i) & iSel(:,i),i,1),...
        clusterPos(iSign_neg(:,i) & iSel(:,i),i,2),10,cols(3,:),...
        'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',.5);%,'MarkerEdgeAlpha',.3);
    
    axis square
    axis equal
    axis off
    
    %xlim([0,5500]); ylim([-5000,1000]);
    xlim([2700,5400]); ylim([-3700, -1000]);
    daspect([1 1 1])
    
end

linkaxes(ax,'xy');

%% plot grid of average proportions per location

plotGrid(squeeze(clusterPos(:,1,1)),squeeze(clusterPos(:,1,2)),clusterVs);


%fpos = [50, 100, 1400, 170];
figure('Position',[50, 100 2500 350]);
tiledlayout(1,size(clusterPs,2));

iSel = clusterPs_contrast>-Inf;

for i=1:size(clusterPs,2)
    
    ax(i) = nexttile;    
    hold on
    
    
    
    %plot atlas boundaries
    cellfun(@(x) cellfun(@(x) plot(1000*x(:,2),1000*x(:,1),'color',[.5 .5 .5]),x,'uni',false), ...
        {bas.dorsal_brain_areas.boundaries_stereotax},'uni', false);
    
    axis square
    axis equal
    axis off
    
    xlim([0,5500]);
    ylim([-5000,1000]);
    daspect([1 1 1])
    
end

linkaxes(ax,'xy');

%% find a good example and plot its p-values

%example_cUID = '49e132f1-cf6f-4ba5-bc13-b67d29d2e03a';
%example_cUID = 'f41e19a7-4545-4e76-b60f-aecc4026ecf8';
example_cUID = '542a8647-0072-4848-bdc5-a2f3b6653d3f';
%example_cUID = '010ef4c7-6a6a-4c17-90cf-88a296e4b4e4';
ix = [AllROIs.clusterUID] == example_cUID;

% plot p-value histograms
figure; 
ixs = find(ix); 
for i=1:length(ixs)
    subplot(1,length(ixs),i); 
    histogram(AllROIs(ixs(i)).tstat_pseudo,20,'edgecolor','none'); 
    hold on; 
    xline(AllROIs(ixs(i)).tstat_empirical,'linewidth',2,'color',[1,0,0]); 
    title(sprintf('p=%.3f',AllROIs(ixs(i)).p)); 
    box off; 
    set(gca,'YColor','none'); 
end

%% sort and plot with interactive viewer

%[~,iSort] = sort(clusterPassiveR(:,end),'descend');

valLastDay = clusterVs(:,end);
[~,iSort] = sort(valLastDay,'descend');

% iSort(1) = find(clusterLookup == '49e132f1-cf6f-4ba5-bc13-b67d29d2e03a');
% iSort(2) = find(clusterLookup == '542a8647-0072-4848-bdc5-a2f3b6653d3f');
% iSort(3) = find(clusterLookup == 'df3f92c0-60bb-450d-9302-d1f21cda0f92');
% iSort(4) = find(clusterLookup == '0527a17f-2b6f-4671-afce-e27c0200bb0f');
% iSort(5) = find(clusterLookup == '55fd2320-d819-4a99-9993-1971714feb49');
% iSort(6) = find(clusterLookup == 'c613d3ea-4c1e-45c4-9af3-27aaaf3db666'); %visually tuned?
% iSort(7) = find(clusterLookup == '1fbdba57-9d29-4f11-abf2-28487c213f24'); %visually tuned and prior tuned?



%valMean = mean(clusterVs,2);
%[~,iSort] = sort(valMean,'descend');

%iSort = 1:size(clusterPETHs,1); %no sorting at all

clusterPETHs_sorted = clusterPETHs(iSort(1:40),:,:,:,:);
cUIDs = clusterLookup(iSort(1:40));

interactiveClusterViewer(clusterPETHs_sorted, T, cUIDs);