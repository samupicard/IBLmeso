% sPaths = IBL_listSessionPaths('root','Y:\Subjects',...
%     'protocol',{'biasedChoiceWorld'}, ...
%     'perfOnEasy',0.9,'trials',490,'mpci',true...
%     );
% sPaths = fliplr(sPaths);

load('canonicalSessions.mat');

delete(gcp('nocreate'));
parpool('Processes');

for i=1:length(sPaths)

    fprintf('\nSession %d/%d: %s\n',i,length(sPaths),sPaths{i});

    try

        %% load neural data, trials table and wheel data

        datpath = sPaths{i};
        splitPath = split(datpath,'\');
        subj = splitPath{end-2};
        date = splitPath{end-1};
        sess = splitPath{end};

        fprintf('\nLoading 2PI traces...');
        Fall = IBL_loadMesoData(subj, date, sess, ...
            'trace', 'spks', 'fast', true);
        badframes = find(Fall.frameQC~=0); % for now, assume all nonzero frameQC is bad
        if isempty(Fall.tr) || isempty(Fall.time)
            warning('Incomplete imaging dataset. Skipping this session.');
            continue
        end
        if size(Fall.tr, 2) ~= numel(Fall.time)
            warning( ...
                ['Unequal number of frames in ROI activity (%d) ' ...
                'and imaging times (%d). Skipping this session.'], ...
                size(Fall.tr, 2), ...
                numel(Fall.time));
            continue
        end
        if ~isempty(badframes)
            warning('non-zero frameQC found - skipping this session.')
            continue;
        end
        %TODO cleanly take out bad frames instead of skipping altogether

        fprintf('Getting event timings... ');
        trialsT = IBL_loadTrialsTable(datpath,'sync','timeline');
        if isempty(trialsT)
            fprintf('no extracted trials found, skipping this session!\n');
            PETH_struct = [];
            return
        end
        fprintf('Done!\n');

        wheelPosition = readNPY(fullfile(datpath,'alf','task_00','_ibl_wheel.position.npy'));
        wheelTimestamps = readNPY(fullfile(datpath,'alf','task_00','_ibl_wheel.timestamps.npy'));

        %select task epoch
        taskStartT = trialsT{1,'intervals_0'}-1;
        taskEndT = trialsT{end,'intervals_1'}+2;
        taskFrames = Fall.time>taskStartT & Fall.time<taskEndT;

        %select good cells
        nROIs = length(Fall.iscell);
        goodCells = Fall.iscell;
        %goodCells = false(size(Fall.iscell)); goodCells(randperm(length(Fall.iscell),50)) = true; %for testing

        Ffull = Fall.tr(goodCells,taskFrames)';
        frameTimes = Fall.time(taskFrames);

        idxFull = Fall.idx(goodCells);
        fovFull = Fall.fov(goodCells);

        %initialize
        deUniqueFull = nan(length(Fall.iscell),1);

        %% loop over FOVs

        fovIxs = unique(fovFull);
        nFOVs = length(fovIxs);

        for iFOV = 1:nFOVs

            fovIx = fovIxs(iFOV);
            fovName = sprintf('FOV_0%d',fovIxs(iFOV));
            F = Ffull(:,fovFull==fovIx);
            idx = idxFull(fovFull==fovIx);

            fprintf('\nFitting %s:\n',fovName);

            %% fit reduced rank models to get kernel selectivity
            rankR = 10;

            results = fitIBLKernelSelectivity( ...
                F, ...
                frameTimes, ...
                trialsT, ...
                wheelPosition, ...
                wheelTimestamps, ...
                rankR, ...
                contrastExponent=0.5, ...
                wheelSpeedSmoothing=0.1, ...
                inTrialMask=true, ...
                nFolds=3, ...
                selectivityThreshold=0.01);

            nPreds  = length(results.predictorNames);
            deUnique = nan(sum(Fall.fov==fovIx),nPreds); %full size (all ROIs)
            deFull = nan(sum(Fall.fov==fovIx),1);

            deUnique(idx+1,:) = results.deUnique;
            deUnique(deUnique<-0.5) = -.5; %HACK to get nothing below -0.5

            deFull(idx+1) = results.full.devianceExplained; % Store the deFull values for the current FOV
            deFull(deFull<-0.5) = -.5; %HACK to get nothing below -0.5

            fprintf('Saving outputs..')
            deFullName = 'mpciROIs.deFull.npy';
            writeNPY(deFull,fullfile(datpath,'alf',fovName,deFullName));

            for iP = 1:nPreds
                deUniqueName = sprintf('mpciROIs.deUnique_%s.npy',results.predictorNames(iP));
                writeNPY(deUnique(:,iP),fullfile(datpath,'alf',fovName,deUniqueName));

            end

            fprintf('. Done!\n');

        end

    catch
       fprintf('\nError, skipping!\n');
    end

end