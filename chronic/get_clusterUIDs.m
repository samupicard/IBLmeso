function clusterUIDs = get_clusterUIDs(subject,date,session)

%loads clusterUIDs from all FOVs of a given session path, and concatenates
%them into a big array of strings.
%if clusterUIDs don't exist for a given FOV, fills with empty strings
%instead
%
%written by Samuel Picard

if nargin==1 %assume we provided a session path
    sessionpath = subject;
    if strcmp(sessionpath(end),'\')
        sessionpath = sessionpath(1:end-1);
    end
    splitPath = split(sessionpath,'\');
    subject = splitPath{end-2};
    date = splitPath{end-1};
    session = splitPath{end};
end

root = 'Y:\Subjects\';
datpath = fullfile(root,subject,date,session);

%fns = dir(fullfile(datpath,'alf','FOV*','mpciROIs.clusterUIDs.csv')); %needs to contain clusterUIDs file
fns = dir(fullfile(datpath,'alf','FOV*','mpciROIs.mpciROITypes.npy')); %needs to contain suite2p outputs

clusterUIDs = [];

for iFOV = 1:length(fns)
    
    reg_data_path = fns(iFOV).folder;
    
    if isfile(fullfile(fns(iFOV).folder,'mpciROIs.clusterUIDs.csv'))
        % Read the CSV file (assumes one column, no header)
        cluster_ids_1FOV = read_uid_csv(fullfile(fns(iFOV).folder,'mpciROIs.clusterUIDs.csv'));
    else
        % make array of empty strings
        roitypes = readNPY(fullfile(fns(iFOV).folder,'mpciROIs.mpciROITypes.npy'));
        cluster_ids_1FOV = strings(length(roitypes),1);
    end
    clusterUIDs = [clusterUIDs; cluster_ids_1FOV];
    
end