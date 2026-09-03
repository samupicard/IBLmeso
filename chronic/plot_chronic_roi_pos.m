function XY_all = plot_chronic_roi_pos(sessionPath, chronicUIDs, opt)
% plot_chronic_roi_pos
%
% Inputs:
%   sessionPath  - path to session root OR session/alf
%   chronicUIDs  - string / cellstr / char array of UIDs to plot
%
% Output:
%   XY_all       - concatenated [N x 2] positions plotted

if nargin < 3, opt = struct(); end

% optional region filter inputs
if ~isfield(opt,'region'), opt.region =  {'VIS*','RSP*','SS*','AUD*','MO*'}; end % acronyms/tokens
if ~isfield(opt,'st'),     opt.st = ''; end               % pre-loaded structure tree table
if ~isfield(opt,'stPath'), opt.stPath = "C:\Users\Samuel\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv"; end           % optional if you want to load by path


% Normalize UID list
chronicUIDs = string(chronicUIDs(:));

% Resolve alf path
sessionPath = string(sessionPath);
if endsWith(sessionPath, filesep + "alf") || ...
   endsWith(sessionPath, "/alf") || ...
   endsWith(sessionPath, "\alf")
    alfPath = sessionPath;
else
    alfPath = fullfile(sessionPath, 'alf');
end

if ~isfolder(alfPath)
    error('alf folder not found at: %s', alfPath);
end

% Find FOVs
FOV_dirs = dir(fullfile(alfPath, 'FOV*'));
FOV_dirs = FOV_dirs([FOV_dirs.isdir]);

if isempty(FOV_dirs)
    error('No FOV directories found in %s', alfPath);
end

XY_all = [];
loc_all = [];

for i = 1:numel(FOV_dirs)

    fovPath = fullfile(FOV_dirs(i).folder, FOV_dirs(i).name);
    
    locFile = fullfile(fovPath, 'mpciROIs.brainLocationIds_ccf_2017_estimate.npy');
    posFile = fullfile(fovPath, 'mpciROIs.mlapdv_estimate.npy');
    uidFile = fullfile(fovPath, 'mpciROIs.clusterUIDs.csv');

    if ~isfile(posFile) || ~isfile(uidFile)
        continue
    end
    
    roi_loc = readNPY(locFile);                 % nROIs x 1
    roi_pos = readNPY(posFile);                 % nROIs x 3
    roi_cUID = string(read_uid_csv(uidFile));   % nROIs x 1
    roi_cUID = roi_cUID(:);

    if size(roi_pos,2) < 2
        warning('roi_pos has <2 columns in %s; skipping.', fovPath);
        continue
    end

    if size(roi_pos,1) ~= numel(roi_cUID)
        warning('UID/pos size mismatch in %s; skipping.', fovPath);
        continue
    end

    % Keep only requested chronic UIDs
    keep = roi_cUID ~= "" & ismember(roi_cUID, chronicUIDs);

    XY_all = [XY_all; double(roi_pos(keep,1:2))]; %#ok<AGROW>
    loc_all = [loc_all; double(roi_loc(keep))];
end

% Get structure tree for label mapping
if isfield(opt,'st') && ~isempty(opt.st)
    st = opt.st;
elseif isfield(opt,'stPath') && strlength(string(opt.stPath)) > 0
    st = loadStructureTree(opt.stPath);
else
    error('Structure tree required for region acronym labeling.');
end
rid_all = double(st.id);
acr_all = string(st.acronym);

% Map region IDs -> ix and acronyms
regionList = unique(loc_all);
regionLabels = strings(size(regionList));
region_ix = nan(1,numel(regionList));
for i = 1:numel(regionList)
    idx = find(rid_all == regionList(i), 1);
    if ~isempty(idx)
        regionLabels(i) = acr_all(idx-1); %hack to remove '1'
    else
        regionLabels(i) = string(regionList(i)); % fallback
    end
    region_ix(i) = idx;
end


% Plot
figure('Color','w','Position',[200,400,280,210]);

hold on;

% Plot points with fixed colors
[tfRid, idxRow] = ismember(double(loc_all), rid_all); 
cmap = turbo(350);
scatter(XY_all(:,1), XY_all(:,2), 5, cmap(idxRow-20, :));

%plot atlas boundaries on top
bas = aratopdown.atlas.build_topdown;
cellfun(@(x) cellfun(@(x) plot(1000*x(:,2),1000*x(:,1),'color',[.5 .5 .5]),x,'uni',false), ...
    {bas.dorsal_brain_areas.boundaries_stereotax},'uni', false);

axis equal;
xlim(mean(XY_all(:,1))+[-1700,1700]); xlabel('ML')
ylim(mean(XY_all(:,2))+[-1700,1700]); ylabel('AP')
daspect([1 1 1])
            
xlabel('ML');
ylabel('AP');

% Title
title(sprintf('Chronic ROIs (N=%d)', size(XY_all,1)));

end