function Fall = IBL_loadMesoData(subject,date,session,varargin)
%IBL_loadMesoData loads IBL alf data of a given session, and puts them into
%a single Fall struct (analogous to suite2p). Loads all FOVs or a selection.
%
% NOTE (speed refactor):
% - Avoids repeated concatenation inside loop (uses cell buffers + vertcat once)
% - Vectorizes ROI->pixel mapping and timeshift extraction
% - Uses ismember for CCF id -> structure index mapping
% - Reduces filesystem calls where possible
%
% readNPY returns arrays as (time x roi) for these files.

if nargin==1 %assume we provided a session path
    sessionpath = subject;
    splitPath = split(sessionpath,'\');
    subject = splitPath{end-2};
    date = splitPath{end-1};
    session = splitPath{end};
end

p = inputParser;
addParameter(p, 'fov','',@(x) ischar(x) || isnumeric(x));
addParameter(p, 'trace', '', @(x) ischar(x) || (isnumeric(x) && isempty(x)));
addParameter(p, 'fast', false, @islogical);
parse(p,varargin{:});
[fov, trace, fast] = deal(p.Results.fov, p.Results.trace, p.Results.fast);

if ischar(fov)
    fovID = fov;
else
    fovID = sprintf('FOV_0%d',fov-1);
end

cellClassifier_thresh = .5;

root = 'Y:\Subjects\';
datpath = fullfile(root,subject,date,session);

% load Allen CCF structure table
if ~fast
    locroot = 'C:\Users\Samuel\Documents\';
    CCFdir = fullfile(locroot,'GitHub','allenCCF');
    st = loadStructureTree(fullfile(CCFdir,'structure_tree_safe_2017.csv'));
    % cache IDs for fast mapping
    stIds = [st.id];
end

%% find FOV folders (by presence of mpci.ROIActivityF.npy)
if isempty(fov)
    fns = dir(fullfile(datpath,'alf','FOV*','mpci.ROIActivityF.npy'));
else
    fns = dir(fullfile(datpath,'alf',fovID,'mpci.ROIActivityF.npy'));
end

% metadata (only needed for some time/pos fallbacks)
meta = [];
if ~fast
    tiff_fn = dir(fullfile(datpath,'raw_imaging_*'));
    try
        meta_fn = dir(fullfile(tiff_fn(1).folder,tiff_fn(1).name,'*meta_NEW.mat'));
        if isempty(meta_fn)
            meta_fn = dir(fullfile(tiff_fn(1).folder,tiff_fn(1).name,'*_2P_*.mat'));
            if size(meta_fn,1)>1
                meta_fn = dir(fullfile(tiff_fn(1).folder,tiff_fn(1).name,'*_2P_*tif.mat'));
            end
            load(fullfile(meta_fn.folder,meta_fn.name),'meta');
        else
            meta = load(fullfile(meta_fn.folder,meta_fn.name));
        end
        % alternative json
        jsontxt = fileread(fullfile(tiff_fn(1).folder,tiff_fn(1).name,'_ibl_rawImagingData.meta.json'));
        meta = jsondecode(jsontxt);
    catch
        warning('could not find valid meta-data file!')
    end
end

%% preallocate per-FOV buffers to avoid repeated concatenation
nFOV = numel(fns);
if nFOV==0
    Fall = struct();
    warning('No FOVs found.');
    return
end

maxnr = 11;
if nFOV > maxnr
    warning(sprintf('more than %d FOVs found, skipping this session for now',maxnr))
    Fall = struct();
    return
end

wantSingleTrace = ~isempty(trace) && ischar(trace);

% trace buffers
F_c    = cell(nFOV,1);
Fneu_c = cell(nFOV,1);
spks_c = cell(nFOV,1);
tr_c   = cell(nFOV,1);

% meta buffers
iscell_c  = cell(nFOV,1);
fov_c     = cell(nFOV,1);
idx_c     = cell(nFOV,1);
cluster_c = cell(nFOV,1);

% optional buffers
pos_c    = cell(nFOV,1);
annot_c  = cell(nFOV,1);
ccfid_c  = cell(nFOV,1);
tshift_c = cell(nFOV,1);

times = [];
times0 = [];
frameQC = [];

error_flag = false(nFOV,1);

fprintf('FOV ');
nCharsPrinted = 0;

for iFOV = 1:nFOV
    fprintf(repmat('\b', 1, nCharsPrinted))
    nCharsPrinted = fprintf('%d/%d ..', iFOV, nFOV);

    reg_data_path = fns(iFOV).folder;
    planeNm = sprintf('plane%s',reg_data_path(end));
    %reg_data_path_masknmf = fullfile(datpath,'suite2p',planeNm,'masknmf_output');
    reg_data_path_masknmf = reg_data_path;
    try
        %% load traces
        if strcmp(trace,'spks')
            % mpci.ROIActivityDeconvolved.npy is (time x roi)
            tr = readNPY(fullfile(reg_data_path_masknmf,'mpci.ROIActivityDeconvolved.npy'))';
        elseif strcmp(trace,'dFF')
            F    = readNPY(fullfile(reg_data_path,'mpci.ROIActivityF.npy'))';
            Fneu = readNPY(fullfile(reg_data_path,'mpci.ROINeuropilActivityF.npy'))';
            tr = F - 0.7*Fneu;
        else
            F    = readNPY(fullfile(reg_data_path,'mpci.ROIActivityF.npy'))';
            Fneu = readNPY(fullfile(reg_data_path,'mpci.ROINeuropilActivityF.npy'))';
            tr   = readNPY(fullfile(reg_data_path,'mpci.ROIActivityDeconvolved.npy'))';
        end

        cellClassifier = readNPY(fullfile(reg_data_path,'mpciROIs.cellClassifier.npy'));
        iscell = cellClassifier > cellClassifier_thresh;
        %iscell = true(size(tr,1),1);

        nROI = numel(iscell);

        %% times
        times = readNPY(fullfile(reg_data_path,'mpci.times.npy'));
        if isempty(times0)
            times0 = times;
        end
        FOV_timeshift = mean(times-times0);
        
        %% frameQC (only once, from first valid FOV)
        if isempty(frameQC)

            % quick hack to deal with time vector being neuralFrames (instead of volumeFrames)
            try
                if size(tr,2) < size(times,1) && isfield(meta,'scanImageParams') && ...
                        isfield(meta.scanImageParams,'hStackManager') && ...
                        isfield(meta.scanImageParams.hStackManager,'zs') && ...
                        size(meta.scanImageParams.hStackManager.zs,1) > 1
                    zstep = size(meta.scanImageParams.hStackManager.zs,1);
                    times = times(1:zstep:end);
                end
            catch
                % ignore
            end

            qc1 = fullfile(reg_data_path,'mpci.mpciFrameQC.npy');
            qc2 = fullfile(reg_data_path,'mpci.badFrames.npy');
            if isfile(qc1)
                frameQC = readNPY(qc1);
            elseif isfile(qc2)
                frameQC = uint8(readNPY(qc2));
            else
                frameQC = zeros(size(tr,2),1,'uint8');
            end
        end

        %% optional: positions / annotation
        stackPos = readNPY(fullfile(reg_data_path_masknmf,'mpciROIs.stackPos.npy')); % (nROI x 2) [x y]?
        if ~fast
            % try preferred precomputed files first (avoid expensive computation)
            roiCCFName   = dir(fullfile(reg_data_path,'mpciROIs.brainLocationIds_ccf_2017_estimate*'));
            roiMLAPDVName = dir(fullfile(reg_data_path,'mpciROIs.mlapdv_estimate*'));

            if ~isempty(roiCCFName) && ~isempty(roiMLAPDVName)
                ccf_id = readNPY(fullfile(reg_data_path,roiCCFName(1).name));

                % map ccf_id -> row index in structure tree (fast)
                [tf, loc] = ismember(ccf_id, stIds);
                annot = nan(size(ccf_id));
                annot(tf) = loc(tf);

                pos = readNPY(fullfile(reg_data_path,roiMLAPDVName(1).name));

            else
                roiAnnotName = dir(fullfile(reg_data_path,'mpciROIs.brainLocationIds*'));
                if ~isempty(roiAnnotName) && ~isempty(roiMLAPDVName)
                    pos = readNPY(fullfile(reg_data_path,roiMLAPDVName(1).name));
                    annot = readNPY(fullfile(reg_data_path,roiAnnotName(1).name))';
                    try
                        ccf_id = st{annot,'id'};
                    catch
                        ccf_id = nan(size(annot));
                    end
                else
                    % compute from pixel maps if present (vectorized)
                    pixelAnnotName = dir(fullfile(reg_data_path,'mpciMeanImage.brainLocationIds*'));
                    pixelMLAPDVName = dir(fullfile(reg_data_path,'mpciMeanImage.mlapdv*'));

                    pixelMLAPDV = [];
                    pixelAnnot  = [];

                    if ~isempty(pixelAnnotName) && ~isempty(pixelMLAPDVName)
                        pixelMLAPDV = readNPY(fullfile(reg_data_path,pixelMLAPDVName(1).name));

                        if ~isempty(regexpi(pixelAnnotName(1).name,'ccf_2017','once'))
                            pixelCCF = readNPY(fullfile(reg_data_path,pixelAnnotName(1).name));
                            % map each pixel CCF id to structure index
                            [tfPix, locPix] = ismember(pixelCCF(:), stIds);
                            pixelAnnot = nan(size(pixelCCF),'like',double(1));
                            tmp = nan(numel(pixelCCF),1);
                            tmp(tfPix) = locPix(tfPix);
                            pixelAnnot(:) = tmp;
                        else
                            pixelAnnot = readNPY(fullfile(reg_data_path,pixelAnnotName(1).name));
                        end
                    else
                        % fallback to meta if available
                        try
                            pixelMLAPDV = meta.FOV(iFOV).pixelMLAPDV;
                            pixelAnnot  = meta.FOV(iFOV).pixelAnnot;
                        catch
                            pixelMLAPDV = [];
                            pixelAnnot  = [];
                        end
                    end

                    n = size(stackPos,1);
                    pos   = nan(n,3);
                    annot = nan(n,1);

                    if ~isempty(pixelAnnot) && ~isempty(pixelMLAPDV)
                        % stackPos(:,1)=x, stackPos(:,2)=y, but indexing is (y,x)
                        x = stackPos(:,1);
                        y = stackPos(:,2);

                        % bounds-safe mask
                        H = size(pixelAnnot,1);
                        W = size(pixelAnnot,2);
                        inb = x>=1 & x<=W & y>=1 & y<=H;

                        lin = nan(n,1);
                        lin(inb) = sub2ind([H W], y(inb), x(inb));

                        annot(inb) = pixelAnnot(lin(inb));

                        % pixelMLAPDV is H x W x 3
                        HW = H*W;
                        lin2 = lin(inb);
                        pos(inb,1) = pixelMLAPDV(lin2 + 0*HW);
                        pos(inb,2) = pixelMLAPDV(lin2 + 1*HW);
                        pos(inb,3) = pixelMLAPDV(lin2 + 2*HW);
                    end

                    try
                        ccf_id = st{annot,'id'};
                    catch
                        ccf_id = nan(size(annot));
                    end
                end
            end
        end

        %% optional: timeshift (vectorized)
        %if ~fast
            tsFile = fullfile(reg_data_path,'mpciStack.timeshift.npy');
            if isfile(tsFile)
                timeshift_per_line = readNPY(tsFile)'; % row/col doesn't matter for linear indexing
            else
                % size based on max x in stackPos (line index)
                if exist('stackPos','var') && ~isempty(stackPos)
                    L = max(512, max(stackPos(:,1)));
                else
                    L = 512;
                end
                timeshift_per_line = zeros(L,1);
            end

            if exist('stackPos','var') && ~isempty(stackPos)
                xline = round(stackPos(:,1));
                % guard if xline exceeds vector length
                if max(xline) > numel(timeshift_per_line)
                    timeshift_per_line(end+1:max(xline),1) = 0; %#ok<AGROW>
                end
                timeshift = timeshift_per_line(xline) + FOV_timeshift;
                timeshift = timeshift(:);
            else
                timeshift = nan(nROI,1);
            end
        %end

        %% clusterUIDs
        clusterUID_path = fullfile(reg_data_path,'mpciROIs.clusterUIDs.csv');
        if isfile(clusterUID_path)
            cluster_ids = read_uid_csv(clusterUID_path);
        else
            cluster_ids = strings(nROI,1);
        end

        %% store into per-FOV buffers (NO concatenation here)
        if wantSingleTrace
            tr_c{iFOV} = tr;
        else
            F_c{iFOV}    = F;
            Fneu_c{iFOV} = Fneu;
            spks_c{iFOV} = tr;
        end

        iscell_c{iFOV}  = iscell;
        fov_c{iFOV}     = repmat(iFOV-1, size(iscell));
        idx_c{iFOV}     = (0:nROI-1).';
        cluster_c{iFOV} = cluster_ids;

        if ~fast
            pos_c{iFOV}    = pos;
            annot_c{iFOV}  = annot;
            ccfid_c{iFOV}  = ccf_id;
        end
        tshift_c{iFOV} = timeshift;


    catch
        error_flag(iFOV) = true;
        warning(sprintf('FOV_%02d has an irregularity, skipping!',iFOV-1))
        continue
    end
end

% if all FOVs failed, bail
if all(error_flag)
    Fall = struct();
    return
end

% keep only successful FOVs
good = ~error_flag;

if wantSingleTrace
    tr_c = tr_c(good);
else
    F_c    = F_c(good);
    Fneu_c = Fneu_c(good);
    spks_c = spks_c(good);
end

iscell_c  = iscell_c(good);
fov_c     = fov_c(good);
idx_c     = idx_c(good);
cluster_c = cluster_c(good);

if ~fast
    pos_c   = pos_c(good);
    annot_c = annot_c(good);
    ccfid_c = ccfid_c(good);
end
tshift_c= tshift_c(good);

%% build Fall once
if fast
    if wantSingleTrace
        Fall = struct('tr',[],'iscell',[],'time',[],'frameQC',[],'fov',[],'idx',[],'clusterUID',[]);
    else
        Fall = struct('F',[],'Fneu',[],'spks',[],'iscell',[],'time',[],'frameQC',[],'fov',[],'idx',[],'clusterUID',[]);
    end
else
    if wantSingleTrace
        Fall = struct('tr',[],'iscell',[],'time',[],'pos',[],'annot',[],'ccf_2017_id',[],'timeshift',[],'frameQC',[],'fov',[],'idx',[],'clusterUID',[]);
    else
        Fall = struct('F',[],'Fneu',[],'spks',[],'iscell',[],'time',[],'pos',[],'annot',[],'ccf_2017_id',[],'timeshift',[],'frameQC',[],'fov',[],'idx',[],'clusterUID',[]);
    end
end

if wantSingleTrace
    Fall.tr = vertcat(tr_c{:});
else
    Fall.F    = vertcat(F_c{:});
    Fall.Fneu = vertcat(Fneu_c{:});
    Fall.spks = vertcat(spks_c{:});
end

Fall.iscell     = vertcat(iscell_c{:});
Fall.time       = times;
Fall.frameQC    = frameQC;
Fall.fov        = vertcat(fov_c{:});
Fall.idx        = vertcat(idx_c{:});
Fall.clusterUID = vertcat(cluster_c{:});

if ~fast
    Fall.pos         = vertcat(pos_c{:});
    Fall.annot       = vertcat(annot_c{:});
    Fall.ccf_2017_id = vertcat(ccfid_c{:});
end
Fall.timeshift   = vertcat(tshift_c{:});

%% consistency checks / patches (preserve your existing behavior)
if wantSingleTrace
    traceLength = size(Fall.tr,2);
else
    traceLength = size(Fall.F,2);
end

if size(Fall.time,1) < traceLength
    warning('Nr of frames in time vector is smaller than nr of frames in neural data!')
    if strcmp(subject,'SP035') && strcmp(date,'2023-04-21') && strcmp(session,'001')
        fprintf('EXCEPTION: Patching up mismatching frame nrs by chopping off the beginning..')
        if wantSingleTrace
            Fall.tr = Fall.tr(:, end-size(Fall.time,1)+1:end);
        else
            Fall.F    = Fall.F(:,    end-size(Fall.time,1)+1:end);
            Fall.Fneu = Fall.Fneu(:, end-size(Fall.time,1)+1:end);
            Fall.spks = Fall.spks(:, end-size(Fall.time,1)+1:end);
        end
    else
        return
    end
elseif size(Fall.time,1) > traceLength
    warning('Nr of frames in time vector is LARGER than nr of frames in neural data?!')
end

fprintf('. Done!\n');
end
