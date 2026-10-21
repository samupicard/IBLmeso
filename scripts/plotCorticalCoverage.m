% Plot mesoscope coverage and ROI locations.
root = 'Y:\Subjects';
all_subjects = {'SS010'};
test_path = 'SS010\2026-09-24\004'; col = [192 0 0]/255;
protocol = {'trainingChoiceWorld','biasedChoiceWorld'};
chronic_filter = false; %if true, will only plot chronic ROIs once
saveflag = false;

mlapdv_fn = 'mpciMeanImage.mlapdv_estimate.npy';
%load allenCCF labels and colors

CCFdir = 'C:\Users\Samuel\Documents\GitHub\allenCCF\\';
av = readNPY(fullfile(CCFdir,'annotation_volume_10um_by_index.npy'));
st = loadStructureTree(fullfile(CCFdir,'structure_tree_safe_2017.csv'));
load(fullfile(CCFdir,'Browsing Functions','allen_ccf_colormap_2017.mat'),'cmap');
bas = aratopdown.atlas.build_topdown;
colmap = hsv(length(all_subjects));
%% plot FOV positions
if true

    posList_name = {};
    posList_pos = [];
    posList_cnt = 0;
    nFOVs = 0;
    nSess = 0;
    coords_means = [];
    for iSubj = 1:length(all_subjects)

        paths = {[root filesep test_path]};
        splitPaths = split(paths,filesep);
        figure('Position',[250   400   850   800]);
        if size(splitPaths,2)==0

            subject = all_subjects{iSubj};
            dates = [''];
            sessions = [''];
        elseif size(splitPaths,2)==1

            subjects = splitPaths(end-2);
            dates = splitPaths(end-1);
            sessions = splitPaths(end);
        else

            subjects = splitPaths(:,:,end-2);
            dates = splitPaths(:,:,end-1);
            sessions = splitPaths(:,:,end);
        end
        nsessions = 0;
        coords_means = [];
        allPos = double.empty(0,3);
        for iSess = 1:length(sessions)

            subject = subjects{iSess};
            date = dates{iSess};
            session = sessions{iSess};
            datloc = fullfile(subject,date,session);
            datpath = fullfile(root,datloc);

            %OPTION 1: plot FOV positions in 'degrees from centre' using the raw imaging metadata

            if false

                meta_dir = dir(fullfile(datpath,'raw_imaging_data*','*tif.mat'));
                load(fullfile(meta_dir(1).folder, meta_dir(1).name));
                figure;
                hold on;
                cols = jet(size(meta.FOV,2));
                for i=1:size(meta.FOV,2)

                    d = meta.FOV(i).MM;
                    coords = 1000*[d.topLeft; d.topRight; d.bottomRight; d.bottomLeft]';
                    patch('XData',coords(1,:), 'YData', coords(2,:), 'EdgeColor',cols(i,:),'FaceColor',[1, 1, 1], 'FaceAlpha', 0);
                    text(mean(coords(1,:)),mean(coords(2,:)),sprintf('%02d',i-1),'color',cols(i,:),'horizontalalignment','center');
                end
                xlim([-2500,2500]); ylim([-3000,2000]);
                xlabel('ML'); ylabel('AP')

                daspect([1 1 1])

                subplot(2,1,2);
                xline(0,'--')

                yline(0,'--')

                cols = jet(size(meta.FOV,2));
                for i=1:size(meta.FOV,2)

                    d = meta.FOV(i).Deg;
                    coords = [d.topLeft; d.topRight; d.bottomRight; d.bottomLeft]';
                    patch('XData',coords(1,:), 'YData', coords(2,:), 'EdgeColor',cols(i,:),'FaceColor',[1, 1, 1], 'FaceAlpha', 0);
                    text(mean(coords(1,:)),mean(coords(2,:)),sprintf('%02d',i-1),'color',cols(i,:),'horizontalalignment','center');
                end
                set(gca,'YDir','reverse')

                xlim([-15,15])
                ylim([-15,15])

                xlabel('X (deg)')
                ylabel('Y (deg)')

                daspect([1 1 1])

                %OPTION 2: plot FOV positions using the mpci data

            else

                nexttile;
                hold on;
                FOV_dirs = dir(fullfile(datpath,'alf','FOV*'));
                errorflag = false(1,size(FOV_dirs,1));
                cols = jet(size(FOV_dirs,1));
                for i=1:size(FOV_dirs,1)

                    mlapdv_pos = readNPY(fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, mlapdv_fn));
                    patch_col = [192 0 0]/255;
                    if ndims(mlapdv_pos)==2

                        mlapdv_pos = permute(reshape(mlapdv_pos, 512, 512, 3), [2 1 3]);
                    end
                    [V, F] = get_valid_corners(mlapdv_pos);
                    coords = V(:,1:2)';
                    coords_means = [coords_means mean(coords,2)];
                    %OPTION 1

                    meanImg = readNPY(fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, 'mpciMeanImage.images.npy'));
                    valid = all(isfinite(mlapdv_pos), 3);
                    X = mlapdv_pos(:,:,1);
                    Y = mlapdv_pos(:,:,2);
                    Z = mlapdv_pos(:,:,3);
                    Ivalid = meanImg(valid);
                    Ivalid_norm = (Ivalid - min(Ivalid)) ./ (max(Ivalid) - min(Ivalid));
                    scatter(X(valid), Y(valid), 5, meanImg(valid), '.'); colormap gray; caxis(prctile(meanImg(valid), [5 99]));
                    allPos = cat(1,allPos,reshape(mlapdv_pos,512*512,3));
                    %OPTION 2

                    text(mean(coords(1,:)),mean(coords(2,:)),sprintf('%02d',i-1),'color',cols(iSess,:),'horizontalalignment','center');
                    posList_cnt = posList_cnt + 1;
                    posList_name{posList_cnt} = sprintf('%s/FOV%02d',datloc,i-1);
                    posList_pos = [posList_pos, coords(:,1)];
                    errorflag(i) = true;
                    nFOVs = nFOVs + 1;
                end
                if any(errorflag)

                    nsessions = nsessions+1;
                    nSess = nSess+1;
                end
            end
            cellfun(@(x) cellfun(@(x) plot(1000*x(:,2),1000*x(:,1),'color',[.5 .5 .5]),x,'uni',false), ...
                {bas.dorsal_brain_areas.boundaries_stereotax},'uni', false);

            xlim([0 5200]); ylim([-5100 3600]);
            xlabel('ML'); ylabel('AP');
            daspect([1 1 1])

            set(gca,'xtick',[],'ytick',[]);
            axis square off

            title(mlapdv_fn(15:end-4),'Interpreter','none')

        end
        set(gca,'xtick',[],'ytick',[]);
        axis square off

        daspect([1 1 1])

    end
end

%% plot position of all neurons (colored by area)
%load Allen CCF structure table

CCFdir = 'C:\Users\Samuel\Documents\GitHub\allenCCF\\';
av = readNPY(fullfile(CCFdir,'annotation_volume_10um_by_index.npy'));
st = loadStructureTree(fullfile(CCFdir,'structure_tree_safe_2017.csv'));
load(fullfile(CCFdir,'Browsing Functions','allen_ccf_colormap_2017.mat'),'cmap');
cmap_ccf = cmap;
cmap_rand = [hsv(370); .5, .5 ,.5]; %this should cover all dorsal-posterior cortical areas

%keep info about coverage

dorsal_cortex_regions = {'VISp','VISpm','VISam','VISa','VISrl','VISal','VISl','RSPagl','RSPd','SSp-bfd','SSp-tr','SSp-ul','SSp-un','SSp-ll','SSp-n','SSs','MOp','AUDd','MOs','SSp-m','AUDp','AUDpo','TEa','VISli','VISpor','VISpl'}; %only optically accessible regions

coverage_by_sessions = zeros(1,length(dorsal_cortex_regions));
coverage_by_animals = zeros(1,length(dorsal_cortex_regions));
coveragePairs_by_sessions = zeros(length(dorsal_cortex_regions));
coverage_by_neurons = zeros(1,length(dorsal_cortex_regions));
nNeurons_tot = 0;
nROIs_tot = 0;
nSess_tot = 0;
figure('Name','AllCells','Position',[608   458   657   668]); %[800   200   1100   1000]);
for iSubj = 1:length(all_subjects)

    if chronic_filter

        paths = IBL_listSessionPaths('root',[root filesep all_subjects{iSubj}],'protocol',protocol,'chronic',true,'verbose',true);
        paths = paths(1:19); %for SP072

    else

        paths = IBL_listSessionPaths('root',[root filesep all_subjects{iSubj}],'protocol',protocol,'mpci',true,'verbose',false);
    end
    splitPaths = split(paths,filesep);
    subplot(4,4,iSubj)

    hold on;
    if size(splitPaths,2)==0

        subject = all_subjects{iSubj};
        dates = [''];
        sessions = [''];
    elseif size(splitPaths,2)==1

        subjects = splitPaths(end-2);
        dates = splitPaths(end-1);
        sessions = splitPaths(end);
    else

        subjects = splitPaths(:,:,end-2);
        dates = splitPaths(:,:,end-1);
        sessions = splitPaths(:,:,end);
    end
    %plot atlas boundaries on bottom

    cellfun(@(x) cellfun(@(x) plot(1000*x(:,2),1000*x(:,1),'color',[.5 .5 .5]),x,'uni',false), ...
        {bas.dorsal_brain_areas.boundaries_stereotax},'uni', false);

    roi_pos_all = int64.empty(0,2);
    roi_loc_all = int64.empty(0,1);
    cUIDs_all = string.empty(0,1);
    cUIDs_full = string.empty(0,1);
    roitype_full = false(0,0);
    for iSess = fliplr(1:length(sessions))

        subject = subjects{iSess};
        date = dates{iSess};
        session = sessions{iSess};
        datloc = fullfile(subject,date,session);
        datpath = fullfile(root,datloc);
        FOV_dirs = dir(fullfile(datpath,'alf','FOV*'));
        roi_locs = int64.empty(0,1);
        for i=1:size(FOV_dirs,1)

            try

                roi_pos = readNPY(fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, 'mpciROIs.mlapdv_estimate.npy'));
                roi_loc = readNPY(fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy'));
                roi_type = readNPY(fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, 'mpciROIs.mpciROITypes.npy'));
                cUID_file = fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, 'mpciROIs.clusterUIDs.csv');
                if chronic_filter

                    if isfile(cUID_file)

                        cUIDs = read_uid_csv(cUID_file);
                        cUIDs_empty = cUIDs == "";
                        cUIDs_old = ismember(cUIDs,cUIDs_full) & ~cUIDs_empty; %boolean mask of clusters we've already plotted

                        cUIDs_new = ~ismember(cUIDs,cUIDs_full) & ~cUIDs_empty;
                        mask = cUIDs_new | (cUIDs_empty & logical(roi_type));
                        cUIDs_full = [cUIDs_full; cUIDs(~cUIDs_empty)];
                        cUIDs_all = [cUIDs_all; cUIDs(mask)];
                        roitype_full = [roitype_full; logical(roi_type(~cUIDs_empty))];
                    else

                        warning(sprintf('%s, %s has no cUID file - considering all ROIs classified as cells.',datloc,FOV_dirs(i).name))

                        mask = logical(roi_type);
                        cUIDs_all = [cUIDs_all; strings(sum(mask),1)];
                    end
                else

                    mask = logical(roi_type);
                    cUIDs_all = [cUIDs_all; strings(size(roi_type(mask)))];
                end
                roi_locs = [roi_locs; roi_loc(mask)];
                roi_pos_all = [roi_pos_all;roi_pos(mask,1:2)];
                roi_loc_all = [roi_loc_all;roi_loc(mask)];
                nROIs_tot = nROIs_tot + length(roi_loc);
                nNeurons_tot = nNeurons_tot + sum(mask);
                if i == 1

                    nSess_tot = nSess_tot+1;
                end
            end
        end
        uniqueAreas = unique(roi_locs);
        uniqueAreas = uniqueAreas(uniqueAreas>0);
        for iArea = 1:length(uniqueAreas)

            arealabel = st{find(st.id==uniqueAreas(iArea)),'acronym'}{1}(1:end-1);
            ii = strcmp(dorsal_cortex_regions,arealabel);
            roicnt = sum(roi_locs==uniqueAreas(iArea));
            coverage_by_sessions = coverage_by_sessions + ii;
            coverage_by_neurons = coverage_by_neurons + roicnt*ii;
            for jArea = 1:length(uniqueAreas)

                jj = strcmp(dorsal_cortex_regions,st{find(st.id==uniqueAreas(jArea)),'acronym'}{1}(1:end-1));
                coveragePairs_by_sessions(ii,:) = coveragePairs_by_sessions(ii,:) + jj;
            end
        end
    end
    roi_area_all = cell(size(roi_loc_all));
    roi_col_all = nan(size(roi_loc_all,1),3);
    roi_ndays_all = ones(size(roi_loc_all));
    roi_nCellType_all = nan(size(roi_loc_all));
    cmap_norm = parula(length(sessions));
    for iROI = 1:length(roi_loc_all)

        iC = cUIDs_full==cUIDs_all(iROI);
        nDays = sum(iC); if nDays>=0, roi_ndays_all(iROI)=nDays; end

        if any(iC), roi_nCellType_all(iROI) = sum(roitype_full(iC)); else roi_nCellType_all(iROI)=1; end

        try

            roi_area_all{iROI} = st{find(st.id==roi_loc_all(iROI)),'acronym'}{1};
            roi_col_all(iROI,:) = cmap_ccf(find(st.id==roi_loc_all(iROI)),:);
        end
        roi_ndays_all(iROI) = nDays;
    end
    iSel = roi_nCellType_all>=1 & roi_ndays_all>=0.8*max(roi_ndays_all);
    %plot atlas boundaries on top

    xlim([0,5500]); ylim([-5000,2000]);
    daspect([1 1 1])

    %colorbar stuff

    if chronic_filter

        nTicks = max(roi_ndays_all);
        clim = [0.5,nTicks+0.5];
        colormap(flipud(hot(nTicks)));
        cb = colorbar;
        caxis(clim);
        edges = linspace(clim(1), clim(2), nTicks+1);
        centers = (edges(1:end-1) + edges(2:end)) / 2;
        cb.Ticks = centers;
        cb.TickLabels = string(centers);
        shading flat

        ylabel(cb,'nr of days tracked');
    end
    axis square off

    title([subject ' (n=' num2str(length(roi_loc_all(iSel))) ')'],'Interpreter','none')

    arealabels = {};
    uniqueAreas = unique(roi_loc_all); uniqueAreas = uniqueAreas(uniqueAreas>0);
    for iArea = 1:length(uniqueAreas)

        arealabels{iArea} = st{find(st.id==uniqueAreas(iArea)),'acronym'}{1}(1:end-1);
        coverage_by_animals = coverage_by_animals + strcmp(dorsal_cortex_regions,arealabels{iArea});
    end
end
%% save things
if saveflag
    fn = 'coverage_allCW_eachAnimal_Neurons';
    ff = 'C:\Users\Samuel\Documents\2PI\mesoscope_active\coverage\\';
    %replace impossible pairs of regions with NaNs

    impossible_pairs = getImpossiblePairs;
    for iArea = 1:length(dorsal_cortex_regions)

        regionName = dorsal_cortex_regions{iArea};
        regionName_nodash = erase(regionName,'-');
        for i=1:length(impossible_pairs.(regionName_nodash))

            coveragePairs_by_sessions(iArea,strcmp(dorsal_cortex_regions,impossible_pairs.(regionName_nodash)(i))) = nan;
        end
    end
    %plot region pairs

    load('purples_colormap.mat');  % loads variable `purples`

    gray = [0.75 0.75 0.75];  % light gray for NaNs

    white = [1 1 1];
    new_colormap = [gray; white; purples(round(linspace(20,256,8)),:)];  % prepend gray

    figure;
    colormap(new_colormap);
    data = coveragePairs_by_sessions;
    data(isnan(data)) = -1;
    imagesc(data);
    set(gca,'xtick',1:length(dorsal_cortex_regions),'ytick',1:length(dorsal_cortex_regions),...
        'xticklabels',dorsal_cortex_regions,'yticklabels',dorsal_cortex_regions);

    xtickangle(gca,60);
    caxis([-1.5,8.5]);
    title('Pairs of regions')

    axis image;
    c = colorbar;
    c.Label.String = 'nr of sessions';
    c.Ticks = [-1,0:2:8];
    ticks = c.Ticks;
    tick_labels = cellstr(string(ticks));
    tick_labels{end} = ['≥ ' tick_labels{end}];
    tick_labels{1} = 'impossible';
    c.TickLabels = tick_labels;
    regionName = dorsal_cortex_regions';
    sessionNr = coverage_by_sessions';
    animalNr = coverage_by_animals';
    neuronNr = coverage_by_neurons';
    coverageTable = table(regionName,sessionNr,animalNr,neuronNr);
    writetable(coverageTable,'C:\Users\Samuel\Documents\Python\Data\coverage_allCW\.csv')

end
%% plot position of all neurons in a single figure
%load Allen CCF structure table
if false

    CCFdir = 'C:\Users\Samuel\Documents\GitHub\allenCCF\\';
    av = readNPY(fullfile(CCFdir,'annotation_volume_10um_by_index.npy'));
    st = loadStructureTree(fullfile(CCFdir,'structure_tree_safe_2017.csv'));
    load(fullfile(CCFdir,'Browsing Functions','allen_ccf_colormap_2017.mat'),'cmap');
    cmap_ccf = cmap;
    roi_pos_all = int64.empty(0,2);
    roi_loc_all = int64.empty(0,1);
    figure('Name','AllCells','Position',[257,735,1081,549]);
    hold on;
    for iSubj = 1:length(all_subjects)

        paths = IBL_listSessionPaths('root',[root filesep all_subjects{iSubj}],'protocol',protocol,'raw_imaging',true,'verbose',false);
        splitPaths = split(paths,filesep);
        if size(splitPaths,2)==0

            subject = all_subjects{iSubj};
            dates = [''];
            sessions = [''];
        elseif size(splitPaths,2)==1

            subjects = splitPaths(end-2);
            dates = splitPaths(end-1);
            sessions = splitPaths(end);
        else

            subjects = splitPaths(:,:,end-2);
            dates = splitPaths(:,:,end-1);
            sessions = splitPaths(:,:,end);
        end
        for iSess = 1:length(sessions)

            subject = subjects{iSess};
            date = dates{iSess};
            session = sessions{iSess};
            datloc = fullfile(subject,date,session);
            datpath = fullfile(root,datloc);
            FOV_dirs = dir(fullfile(datpath,'alf','FOV*'));
            for i=1:size(FOV_dirs,1)

                try

                    roi_pos = readNPY(fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, 'mpciROIs.mlapdv_estimate.npy'));
                    roi_loc = readNPY(fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy'));
                    roi_type = readNPY(fullfile(FOV_dirs(i).folder, FOV_dirs(i).name, 'mpciROIs.mpciROITypes.npy'));
                    roi_pos_all = [roi_pos_all;roi_pos(logical(roi_type),1:2)];
                    roi_loc_all = [roi_loc_all;roi_loc(logical(roi_type))];
                    uniqueAreas = unique(roi_loc);
                    for iArea = 1:length(uniqueAreas)

                        arealabel = st{find(st.id==uniqueAreas(iArea)),'acronym'}{1}(1:end-1);
                    end
                end
            end
        end
    end
    roi_area_all = cell(size(roi_loc_all));
    roi_col_all = nan(size(roi_loc_all,1),3);
    for iROI = 1:length(roi_loc_all)

        roi_area_all{iROI} = st{find(st.id==roi_loc_all(iROI)),'acronym'}{1};
        roi_col_all(iROI,:) = cmap_ccf(find(st.id==roi_loc_all(iROI)),:);
    end
    scatter(roi_pos_all(:,1),roi_pos_all(:,2),5,roi_col_all,...
        'Marker','.','MarkerFaceColor','flat','MarkerEdgeColor','flat');%'MarkerFaceAlpha',.3);%,'MarkerEdgeAlpha',.3);

    %plot area labels at centroids

    centroids = [];
    arealabels = {};
    uniqueAreas = unique(roi_loc_all);
    for iArea = 1:length(uniqueAreas)

        centroids(iArea,:) = nanmean(roi_pos_all(roi_loc_all==uniqueAreas(iArea),:),1);
        arealabels{iArea} = st{find(st.id==uniqueAreas(iArea)),'acronym'}{1}(1:end-1);
    end
    if ~isempty(arealabels)

        t = text(centroids(:,1),centroids(:,2),arealabels,'Color',[0 0 0],'FontSize',8,...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

    end
    axis square

    axis equal

    xlabel('ML')

    ylabel('AP')

    xlim([0,5000]);
    ylim([-4500,500]);
    daspect([1 1 1])

    title([num2str(length(all_subjects)),' subjects (n=' num2str(length(roi_loc_all)) ')'],'Interpreter','none')

end