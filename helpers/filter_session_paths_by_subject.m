function sessionPathsOut = filter_session_paths_by_subject(sessionPaths, subjectList)
% Filter session paths by subject IDs.
%
% Inputs
%   sessionPaths : cellstr or string array of full session paths
%   subjectList  : cellstr or string array, e.g. {'SP044','SP052'}
%
% Output
%   sessionPathsOut : filtered sessionPaths

sessionPaths = string(sessionPaths(:));
subjectList  = string(subjectList(:));

keep = false(size(sessionPaths));

for i = 1:numel(sessionPaths)
    sp = sessionPaths(i);

    % Parse subject from .../<subject>/<YYYY-MM-DD>/<sessionID>
    [parentDir, ~] = fileparts(sp);        % remove sessionID
    [subjDir, ~]   = fileparts(parentDir); % remove date
    [~, subjStr]   = fileparts(subjDir);   % get subject

    keep(i) = ismember(string(subjStr), subjectList);
end

sessionPathsOut = sessionPaths(keep);

end