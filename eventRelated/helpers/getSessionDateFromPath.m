function dt = getSessionDateFromPath(datpath)
% getSessionDateFromPath
%
% Extract session date from path like:
%   Y:\Subjects\SP058\2024-07-18\001

dt = NaT;

sp = split(string(datpath), filesep);
for i = 1:numel(sp)
    tok = regexp(char(sp(i)), '^\d{4}-\d{2}-\d{2}$', 'match', 'once');
    if ~isempty(tok)
        dt = datetime(tok, 'InputFormat', 'yyyy-MM-dd');
        return
    end
end
end