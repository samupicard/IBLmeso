root = 'Y:\Subjects';

ROItable = AllROIs;

%prepare some parameters
if isscalar(params.nComp)
    iComp = 1:params.nComp;
else
    iComp = params.nComp;
end
nROIs_all = height(ROItable);   
iROI_good = false(1,nROIs_all);

fR = 5.0761; %this is the majority frame rate at ibl-mesoscope (8 FOVs)
fDur = 1/fR; 

%parse unique sessions from ROItable
[G, subjlist, datelist, sesslist] = findgroups(ROItable.subject, ROItable.date, ROItable.session);
ns = numel(subjlist); %number of unique sessions
cnt = 0;
valid_sess_cnt = 1;
tic

%for each unique session,
for s = 1:ns
    fprintf('Session %d / %d (%.1f%%)', ...
        s, ns, 100*s/ns);
    
    %select the part of the table corresponding to this session
    idx = (G == s);
    rT = ROItable(idx,:);
    nROIs = sum(idx);

    datpath = fullfile(root,string(subjlist(s)),string(datelist(s)),string(sesslist(s)));
    
    %load trials table
    tT = IBL_loadTrialsTable(datpath);
    
    %add a condition called 'stimSide'
    tT.stimSide = sign(tT.contrastDiff);

    %load PETH file
    fn = sprintf('allNeurons.PETHfull_%s.mat',params.evnt);
    load(fullfile(datpath,'alf',fn));
    
    %check that things match
    if nROIs~=size(allPETH_raw_norm,1),
        warning('Nr. of ROIs in ROI table is not the same as in PETH mat file! Skipping this session')
        continue
    end
    
    
    %initialize PETHavg matrix
    if valid_sess_cnt == 1
        trialTypeVals = unique(tT{:,params.trialTypeField});
        trialTypeVals_sel = [trialTypeVals(iComp,:); trialTypeVals(1+length(trialTypeVals)-iComp,:)]';
        T_all = [-fliplr(0:fDur:(-params.twin_all(1)+fDur)),fDur+(0:fDur:params.twin_all(2))]; %needlessly complicated, but this is making sure we have a time vector that runs from just before the window to just after the window
        pval_all = nan(nROIs_all,1);
        tstat_all = nan(nROIs_all,1);
        h_all = false(nROIs_all,1);
        PETHavg_all = nan(nROIs_all,length(trialTypeVals_sel),length(T_all));    
    end
    
    %get all the ROIs that exist in the table
    iROIs = 1:nROIs; %TODO change this to an actual indexing
    
    %compute the mean resps for contrast of interest
    %[meanResps, typeVals, validType] = getMeanResps_allROIs(allPETH_raw_norm(iROIs,:,:), tT, params.trialTypeField, trialTypeVals_sel, {'choice'},5); 
    [meanResps, typeVals, validType] = getMeanResps_allROIs(allPETH_raw_norm(iROIs,:,:), tT, params.trialTypeField, trialTypeVals_sel); %no balancing of sub-conditions
    
    if ~all(validType)
        warning('Not enough trials from each condition-combination. Skipping this session.')
        continue;
    end
    
    iROI_good(idx) = true;

    % --- resample if needed ---
    sameT = numel(T)==numel(T_all) && all(abs(T(:)-T_all(:)) < 0.1); %we tolerate up to 100ms shift
    if ~sameT
        Y = reshape(permute(meanResps, [3 1 2]), numel(T), []);   % nT_sess x (nROIs*nType)
        Yq = interp1(T(:), Y, T_all(:), 'linear', NaN);                    % nT_all x (nROIs*nType)
        meanResps = permute(reshape(Yq, numel(T_all), nROIs, size(meanResps,2)), [2 3 1]);
        fprintf(' [resampled %d->%d]', numel(T), numel(T_all));
    end
    
    %append results
    ix = cnt + iROIs;
    PETHavg_all(ix,:,:) = meanResps;
    pval_all(ix) = rT{iROIs,'p'};
    tstat_all(ix) = rT{iROIs,'tstat_empirical'};
    h_all(ix) = rT{iROIs,'h'};
    
    cnt = cnt + nROIs;
    valid_sess_cnt = valid_sess_cnt+1;
    
    fprintf('\r');
end

%trim unused preallocation
PETHavg_all = PETHavg_all(1:cnt,:,:);
pval_all = pval_all(1:cnt);
tstat_all = tstat_all(1:cnt);
h_all = h_all(1:cnt);

fprintf('\nDone in %.1f minutes.\n', toc/60);
