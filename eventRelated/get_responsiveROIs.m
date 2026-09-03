function [h,pvals] = get_responsiveROIs(datpath, varargin)

%% options
p = inputParser;
p.addParameter('saveToFOV', false, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
saveToFOV = logical(p.Results.saveToFOV);

%TODO separately define the event types, windows etc to use for each satistical test
% params_all = struct( ...
%     'activity_type',   'deconv', ...
%     'trialTypeField',  {'','','choice','','feedbackType','probabilityLeft'}, ...
%     'trialTypeVals',   {[], [], [], [-1,1], [], [-1,1], [0.2,0.8]}, ...
%     'trialTypeFilter', {'', '', '', '', ''}, ...
%     'evnt',            {'stimOn','stimOn','choiceMovement','feedback','stimOn'}, ...
%     'twin_ev',         {[0,0.4],[-0.5,-0.1],[-0.4,0],[-0.4,0],[0,0.4],[-0.5,-0.1]}, ...
%     );

%path logic
splitPath = split(datpath, filesep);
subj = splitPath{end-2};
dt = splitPath{end-1};
sess = splitPath{end};

fprintf('\n%s:\n ',datpath);

%% load task data

fprintf('Getting event timings... ');
trialsT = IBL_loadTrialsTable(datpath,'sync','timeline');
if isempty(trialsT)
    warning('no extracted trials found. Will try to generate trials table from raw bpod events...');
    %try generating the TrialsTable from the bpod time events
    try
        trialsT = IBL_generateTrialsTable(datpath);
    catch
        trialsT = array2table([]);
    end
end
if isempty(trialsT)
    fprintf('no extracted trials found, skipping this session!\n');
    PETH_struct = [];
    return
%elseif size(trialsT,1)<params_all(1).nTrialsMin
%    fprintf('fewer than %d total trials, skipping this session!\n', params_all(1).nTrialsMin);
%    PETH_struct = [];
%    return
end
fprintf('Done!\n');

%% load neural data

%load frame QC, check if there are any non-zero, and throw warning if so
frameQC = readNPY(fullfile(datpath,'alf','FOV_00','mpci.mpciFrameQC.npy'));
if any(frameQC~=0)
    warning('Badframes / non-zero frameQC found. We will ignore trials in those periods.\n');
end

%load data from all FOVs into a single struct (this takes ~30-60s)
fprintf('Loading 2PI traces...');
Fall = IBL_loadMesoData(subj,dt,sess,'trace','spks','fast',true);

if isempty(Fall.tr) || isempty(Fall.time)
    warning('Incomplete datasets, skipping this session!\n');
    PETH_struct = [];
    return
end

if size(Fall.tr,2) ~= size(Fall.time,1)
    warning(sprintf('Unequal nr of frames in mpci.ROIActivity (%d) and in mpci.times (%d). Skipping this session!\n',size(Fall.tr,2),size(Fall.time,1)));
    PETH_struct = [];
    return
end


sig = Fall.tr;              % deconvolved activity
frameTimes = Fall.time;     % relies on existing mpci.times.npy
Fs = 1/median(diff(Fall.time)); % calculate frameRate from time vector
fov = Fall.fov;             % per-ROI fov label (typically 0..N-1)
idx = Fall.idx;             % per-ROI local index within FOV (if provided)

%additionally derive the badframes from frameQC
badframes = find(Fall.frameQC~=0); % for now, assume all nonzero frameQC is bad

%% compute event-locked responses

fprintf('Computing baseline responses: '); 
[Evk_bl, ~, ~, badTrials_bl] = get_evoked_full(sig, frameTimes, find(Fall.frameQC~=0), trialsT, 'stimOn_times', [-0.5,-0.1], 'none');

fprintf('Computing stim-locked responses: '); 
[Evk_stim, ~, ~, badTrials_stim] = get_evoked_full(sig, frameTimes, find(Fall.frameQC~=0), trialsT, 'stimOn_times', [0,0.4], 'none');

fprintf('Computing choiceMovement-locked responses: '); 
[Evk_ch, ~, ~, badTrials_ch] = get_evoked_full(sig, frameTimes, find(Fall.frameQC~=0), trialsT, 'choiceMovement_times', [-0.4,0], 'none');

fprintf('Computing feedback-locked responses: '); 
[Evk_fb, ~, ~, badTrials_fb] = get_evoked_full(sig, frameTimes, find(Fall.frameQC~=0), trialsT, 'feedback_times', [0,0.4], 'none');

badTrials = badTrials_stim | badTrials_bl | badTrials_ch | badTrials_fb;

nTests = 5;

%% compute valid trial indices
trialSel_chType = trialSelection(trialsT,'choice',[-1,1]);
trialSel_chTime = (trialsT.choiceMovement_times - trialsT.stimOn_times)<1.5;
trialSel_ch = any(trialSel_chType,1) & trialSel_chTime' & ~badTrials;

trialSel_fbType = trialSelection(trialsT,'feedbackType',[-1,1]);
trialSel_fbTime = (trialsT.feedback_times - trialsT.stimOn_times)<2;
trialSel_fb = any(trialSel_fbType,1) & trialSel_fbTime' & ~badTrials;

if numel(unique(trialsT.probabilityLeft))==3
    trialSel_blockType = trialSelection(trialsT,'probabilityLeft',[.2 .8]);
    trialSel_block = any(trialSel_blockType,1) & ~badTrials;
    idx_block_t1 = trialSel_block & trialSel_blockType(1,:);
    idx_block_t2 = trialSel_block & trialSel_blockType(2,:);
    nTests = nTests+1;
else %dummy
    trialSel_block = ~badTrials;
    idx_block_t1 = trialSel_block;
    idx_block_t2 = trialSel_block;
end


%% do statistical tests (based on Steinmetz 2019)

stat_names = {
    'signrank_stimOn_vs_baseline'
    'signrank_choiceMovement_vs_baseline_choice'
    'ranksum_choiceMovement_choice'
    'signrank_feedback_vs_baseline_feedbackType'
    'ranksum_feedback_feedbackType'
};

if numel(unique(trialsT.probabilityLeft))==3
    stat_names{end+1} = 'ranksum_baseline_block_probabilityLeft';
end

% sanity check
assert(numel(stat_names) == nTests, 'pval_names length must match nTests');

idx_ch      = trialSel_ch;
idx_ch_t1   = trialSel_ch & trialSel_chType(1,:);
idx_ch_t2   = trialSel_ch & trialSel_chType(2,:);

idx_fb      = trialSel_fb;
idx_fb_t1   = trialSel_fb & trialSel_fbType(1,:);
idx_fb_t2   = trialSel_fb & trialSel_fbType(2,:);

nNeurons = size(Evk_stim,1);
pvals = nan(nNeurons,nTests);
parfor i = 1:nNeurons
    
    pval= nan(1,nTests);
    
    % STIM vs BL
    pval(1) = signrank(Evk_stim(i,:), Evk_bl(i,:), 'method','approximate');

    % CH vs BL (selected trials)
    pval(2) = signrank(Evk_ch(i,idx_ch), Evk_bl(i,idx_ch), 'method','approximate');

    % CH type1 vs type2
    pval(3) = ranksum(Evk_ch(i,idx_ch_t1), Evk_ch(i,idx_ch_t2), 'method','approximate');

    % FB vs BL (selected trials)
    pval(4) = signrank(Evk_fb(i,idx_fb), Evk_bl(i,idx_fb), 'method','approximate');

    % FB type1 vs type2
    pval(5) = ranksum(Evk_fb(i,idx_fb_t1), Evk_fb(i,idx_fb_t2), 'method','approximate');
    
    % BLOCK 0.2 vs 0.8
    if numel(unique(trialsT.probabilityLeft))==3
        pval(6) = ranksum(Evk_bl(i,idx_block_t1), Evk_bl(i,idx_block_t2), 'method','approximate');
    end
    
    pvals(i, :) = pval;
    
end

h = any(pvals<0.05/nTests,2); %bonferroni-corrected alpha value (as in Steinmetz 2019)

%% TODO compute non-parametric d' for choice, feedback, stimSide, block

%d_fb = dprimes_np(Evk_ch(:,trialSel_fb),trialSel_fbType);

%% save results to each respective FOV folder
if saveToFOV
    uFov = unique(fov(:))';

    for k = 1:numel(uFov)
        thisFov = uFov(k);

        % typical IBL: alf/FOV_00, alf/FOV_01, ...
        fovFolder = fullfile(datpath, 'alf', sprintf('FOV_%02d', thisFov));

        if ~exist(fovFolder, 'dir')
            warning('FOV folder not found: %s (skipping)', fovFolder);
            continue
        end

        roiMask = (fov == thisFov);

        if ~any(roiMask)
            warning('No ROIs found for FOV %d (skipping)', thisFov);
            continue
        end

        % subset to ROIs in this FOV
        pvals_fov = pvals(roiMask, :);
        h_fov = h(roiMask);
        
        outName = sprintf('mpciROIs.taskResponsive');
        outPath = fullfile(fovFolder, outName);
        
        T_ps_fov = array2table(pvals_fov, 'VariableNames', stat_names);
            
        try
            writetable(T_ps_fov, [outPath, 'P.tsv'], 'FileType', 'text', 'Delimiter', '\t');
            %writeNPY(pvals_fov, [outPath, 'P.npy']);
            %writeNPY(h_fov, [outPath, '_h.npy']);
        catch ME
            warning('Failed writing %s\n%s', outPath, getReport(ME,'basic'));
        end
        
%         % also drop the names file into the FOV folder for convenience
%         try
%             fid = fopen(fullfile(fovFolder, 'taskResponsiveP.names.tsv'), 'w');
%             for j = 1:numel(stat_names)
%                 fprintf(fid, '%s\n', stat_names{j});
%             end
%             fclose(fid);
%         catch
%             % non-fatal
%         end
    end
end

