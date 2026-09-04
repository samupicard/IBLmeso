function nm = sessNameFromPath(datpath)
% best-effort session label from path ...\subject\date\session
sp = split(string(datpath), filesep);
if numel(sp) >= 3
    nm = sprintf('%s | %s | %s', sp(end-2), sp(end-1), sp(end));
else
    nm = char(datpath);
end
end