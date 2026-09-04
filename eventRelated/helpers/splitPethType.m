function [evnt, subtype] = splitPethType(pethType)
parts = split(string(pethType), "_");
if numel(parts) < 2
    error('Invalid pethType: %s', pethType);
end
evnt = char(parts(1));
subtype = char(strjoin(parts(2:end), "_"));
end