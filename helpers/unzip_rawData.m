root = 'Y:\Subjects\'; %server that has the data
root2 = 'zaru.cortexlab.net\Subjects\';

protocol = 'biasedChoiceWorld';
paths = IBL_listSessionPaths('protocol',protocol,'mpci',true);

%split paths to get subjects, dates, sessions
splitPaths = split(paths,filesep);
subjects = splitPaths(:,:,end-2);
dates = splitPaths(:,:,end-1);
sessions = splitPaths(:,:,end);

for iSess = 1:length(sessions)
    
    subject = subjects{iSess};
    date = dates{iSess};
    session = sessions{iSess};
    
    %get full data path
    datloc = fullfile(subject,date,session);
    datpath = fullfile(root,datloc);
    datpath2 = fullfile(root2,datloc);
    
    fileList = dir(fullfile(datpath,'**','*ROIData.raw.zip'));
    if ~isempty(fileList)
        for iFOV = 1:length(fileList)
            unzip(fullfile(fileList(iFOV).folder,fileList(iFOV).name),fullfile(fileList(iFOV).folder,fileList(iFOV).name(1:end-4)));
        end
    end
    
end