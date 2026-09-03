function paths = IBL_listSessionPaths(varargin)

%IBL_listSessionPaths(Name,Value)
%
%searches all subdirectories of a given path ('root') for sessions containing the
%specified task protocol ('protocol'). Any regular expression can be
%specified as a protocol. Optionally specify whether some datatypes should
%be present.
%
%Inputs:
%'root'            : path in which to search (string. Default: 'Y:\Subjects')
%'protocol'        : protocol(s) to filter (string, or cell of strings). Default: 'biasedChoiceWorld'
%'raw_imaging'     : whether raw_imaging_data folder should be there (boolean). Default: [true,false]
%'alf'             : whether alf folder should be there (boolean). Default: [true,false]
%'mpci'            : whether 'mpci' data should be complete (boolean). Default: [true,false]
%'meta'            : whether metadata file should be there (boolean). Default: [true,false]
%'chronic'         : whether there are chronic clusters (true/false) or min number of tracked cluster (numeric between 0 and 1). Default: [true,false]
%'trials'          : whether there is at least one trials table (boolean). Default: [true,false]
%'fullContrastSet' : whether all expected contrast values must be present. Default: false
%'perfOnEasy'      : minimum fraction correct on trials with abs(contrastDiff)>=0.5
%                    (numeric between 0 and 1, or empty). Default: []
%'verbose'         : whether to print results (boolean). Default: true
%
% NB: if you set 'mpci' to true, you can leave 'raw_imaging', 'alf' and 'meta' to default
%
%Outputs:
%a cell array of full paths to the top folder of the session.

defaultPath = 'Y:\Subjects';
defaultProtocol = '_iblrig_tasks_biasedChoiceWorld';

p = inputParser;
p.addParameter('root', defaultPath, @ischar)
p.addParameter('protocol', defaultProtocol, @(x) ischar(x) || iscell(x))
p.addParameter('raw_imaging', [true,false], @(x) islogical(x) || isnumeric(x));
p.addParameter('mpci', [true,false], @(x) islogical(x) || isnumeric(x));
p.addParameter('alf', [true,false], @(x) islogical(x) || isnumeric(x));
p.addParameter('meta', [true,false], @(x) islogical(x) || isnumeric(x));
p.addParameter('chronic', [true,false], @(x) islogical(x) || (x>=0 && x<=1));
p.addParameter('trials', [true,false], @(x) islogical(x) || isnumeric(x));
p.addParameter('fullContrastSet', false, @islogical)
p.addParameter('perfOnEasy', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=0 && x<=1));
p.addParameter('verbose', true, @islogical)
p.parse(varargin{:})

datpath = p.Results.root;
mpci_flag = logical(p.Results.mpci);
raw_imaging_flag = logical(p.Results.raw_imaging);
alf_flag = logical(p.Results.alf);
meta_flag = logical(p.Results.meta);
chronic_flag = logical(p.Results.chronic);
trials_flag = logical(p.Results.trials);

if ~iscell(p.Results.protocol)
    protocol_template = {p.Results.protocol};
else
    protocol_template = p.Results.protocol;
end

chronic_minFraction = false;
enoughTracked_flag = true;
if any(chronic_flag) && isnumeric(p.Results.chronic)
    chronic_minFraction = p.Results.chronic;
    enoughTracked_flag = false;
end

trials_minN = false;
enoughTrials_flag = true;
if any(trials_flag) && isnumeric(p.Results.trials)
    trials_minN = p.Results.trials;
    enoughTrials_flag = false;
end

fullContrastSet = p.Results.fullContrastSet;
perfOnEasy = p.Results.perfOnEasy;

exclude_folders = {'junk'};

if false % try
    load('Y:\Subjects\allPaths.mat');
    idx = cellfun(@(x) startsWith(x,datpath), d, 'Unif',0);
    d = d([idx{:}]);
else % catch
    d = findFoldersWithRawTaskData(datpath);
end
allSessionDirs = unique(d);

%exclude invalid folders
for i = 1:length(exclude_folders)
    idx = cellfun(@(x) contains(x,exclude_folders{i}), ...
        allSessionDirs, 'Unif',0);
    allSessionDirs = allSessionDirs(~[idx{:}]);
end

iCnt = 0;
paths = {};
protocols = {};
protocols_all = {};
iCnt_chronic = 0;

for iSess = 1:size(allSessionDirs,2)

    rawtasksettings_dir = dir(fullfile(allSessionDirs{iSess}, ...
        'raw_task_data*','_iblrig_taskSettings.raw.json'));
    rawimaging_dir = dir(fullfile(allSessionDirs{iSess}, ...
        'raw_imaging_data*'));
    rawimaging_meta_dir = dir(fullfile(allSessionDirs{iSess}, ...
        'raw_imaging_data*','*meta.json'));
    alf_dir = dir(fullfile(allSessionDirs{iSess},'alf'));
    chronic_dir = dir(fullfile(allSessionDirs{iSess}, ...
        'alf','FOV*','mpciROIs.clusterUIDs.csv'));
    mpciTime_dir = dir(fullfile(allSessionDirs{iSess}, ...
        'alf','FOV*','mpci.times.npy'));
    mpciROI_dir = dir(fullfile(allSessionDirs{iSess}, ...
        'alf','FOV*','mpci.ROI*'));
    trialsT_dir = dir(fullfile(allSessionDirs{iSess}, ...
        'alf','task*','_ibl_trials.table.pqt'));

    % Reset session-specific flags.
    enoughTracked_flag = true;
    enoughTrials_flag = true;
    enoughContrasts_flag = true;
    enoughPerformance_flag = true;

    %check that protocol listed in task settings matches the protocol we are looking for
    p_str_all = {};
    protocol_full = {};

    if ~isempty(rawtasksettings_dir)

        %check that the path contains the files optionally required
        if ismember(~isempty(rawimaging_dir),raw_imaging_flag) ...
                && ismember(~isempty(alf_dir),alf_flag) ...
                && ismember(~isempty(rawimaging_meta_dir),meta_flag) ...
                && ismember(~isempty(mpciROI_dir),mpci_flag) ...
                && ismember(~isempty(mpciTime_dir),mpci_flag) ...
                && ismember(~isempty(chronic_dir),chronic_flag) ...
                && ismember(~isempty(trialsT_dir),trials_flag)

            if chronic_minFraction
                iCnt_chronic = iCnt_chronic+1;
                cUIDs = get_clusterUIDs(allSessionDirs{iSess});

                %index of reference session for fraction tracked
                ref_sess = 1;

                if iCnt_chronic < ref_sess
                    cUIDs_ref = "placeholder";
                elseif iCnt_chronic == ref_sess
                    cellClass_all = [];

                    for iFOV = 1:numel(mpciTime_dir)
                        cellClass_all = [cellClass_all; readNPY(fullfile( ...
                            mpciTime_dir(iFOV).folder, ...
                            'mpciROIs.cellClassifier.npy'))];
                    end

                    cUIDs_ref = cUIDs(cellClass_all>=0.5);
                    cUIDs_ref = cUIDs_ref(cUIDs_ref~="");
                end

                shared_uids = intersect(cUIDs(cUIDs~=""),cUIDs_ref);
                fraction_tracked = numel(shared_uids)/numel(cUIDs_ref);
                enoughTracked_flag = fraction_tracked>=chronic_minFraction;
            end

            if trials_minN
                trials_goCue = readNPY(fullfile( ...
                    trialsT_dir(1).folder, ...
                    '_ibl_trials.goCueTrigger_times.npy'));

                enoughTrials_flag = numel(trials_goCue)>=trials_minN;
            end

            % Both filters use the same trials table and contrastDiff.
            if fullContrastSet || ~isempty(perfOnEasy)
                try
                    trialsT = parquetread(fullfile( ...
                        trialsT_dir(1).folder, ...
                        '_ibl_trials.table.pqt'));

                    contrastRight = trialsT.contrastRight;
                    contrastRight(isnan(contrastRight)) = 0;

                    contrastLeft = trialsT.contrastLeft;
                    contrastLeft(isnan(contrastLeft)) = 0;

                    contrastDiff = contrastRight-contrastLeft;

                    if fullContrastSet
                        requiredContrasts = double([-1,-.25,-.125,-.0625,0,.0625,.125,.25,1]');
                        uniqueContrasts = double(round(unique(contrastDiff)*10000)/10000);
                        uniqueContrasts_no50 = uniqueContrasts(abs(uniqueContrasts)~=0.5); %just in case the 0.5 slipped in
                        enoughContrasts_flag = isequal(uniqueContrasts_no50,requiredContrasts);
                    end

                    if ~isempty(perfOnEasy)
                        easyTrials = abs(contrastDiff)>=0.5;

                        enoughPerformance_flag = any(easyTrials) && ...
                            mean(trialsT.feedbackType(easyTrials)==1) ...
                            >= perfOnEasy;
                    end

                catch
                    enoughContrasts_flag = ~fullContrastSet;
                    enoughPerformance_flag = isempty(perfOnEasy);
                end
            end

            if enoughTracked_flag && enoughTrials_flag && ...
                    enoughContrasts_flag && enoughPerformance_flag

                for iBout = 1:length(rawtasksettings_dir)

                    jsontxt = fileread(fullfile( ...
                        rawtasksettings_dir(iBout).folder, ...
                        rawtasksettings_dir(iBout).name));

                    taskSettings = jsondecode(jsontxt);

                    try
                        protocol = taskSettings.PYBPOD_PROTOCOL;
                        protocol_full = [protocol_full,protocol];
                    catch
                        continue
                    end

                    p_log = [];
                    p_str = {};

                    for i = 1:length(protocol_template)

                        [iFirst,iLast] = regexp( ...
                            protocol,protocol_template{i},'once');

                        taskMatch = ~isempty(iFirst);

                        if taskMatch
                            %check that the end of the strings are really matched
                            if (~strcmp(protocol_template{i}(end),'*') ...
                                    && iLast==length(protocol)) ...
                                    || strcmp(protocol_template{i}(end),'*')

                                p_log = [p_log,true];
                                p_str = [p_str,protocol];
                            end
                        end
                    end

                    if any(p_log)
                        if chronic_minFraction
                            p_str_all = [p_str_all,p_str(p_log), ...
                                sprintf('%.2f tracked',fraction_tracked)];
                        else
                            p_str_all = [p_str_all,p_str(p_log)];
                        end
                    end
                end

                if ~isempty(p_str_all)
                    iCnt = iCnt+1;
                    paths{iCnt} = allSessionDirs{iSess};
                    protocols{iCnt} = strjoin(p_str_all,', ');
                    protocols_all{iCnt} = strjoin(protocol_full,', ');
                end
            end
        end
    end
end

if p.Results.verbose
    fprintf('Found %d valid [',length(paths));

    for i = 1:length(protocol_template)
        fprintf(' %s',protocol_template{i});
    end

    fprintf(' ] sessions:\n');
    fprintf('\n')

    for i = 1:length(paths)
        fprintf('%s\n',paths{i}(end-19:end))
        %fprintf('%s\n',protocols_all{i})
        %fprintf('%s   %s\n',paths{i}(end-19:end),protocols{i})
        %fprintf('%s   %s\n',paths{i},protocols_all{i})
    end
end

end