function [subjects, dates, sessions] = IBL_splitPaths(paths)

%split paths to get subjects, dates, sessions
splitPaths = split(paths,filesep);
if size(splitPaths,2)==1
    subjects = splitPaths(end-2);
    dates = splitPaths(end-1);
    sessions = splitPaths(end);
else
    subjects = splitPaths(:,:,end-2);
    dates = splitPaths(:,:,end-1);
    sessions = splitPaths(:,:,end);
end
