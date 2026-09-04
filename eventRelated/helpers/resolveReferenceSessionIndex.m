function refIdx = resolveReferenceSessionIndex(out, refSpec)
% resolve reference session from numeric index or date string

nSessions = numel(out.S);

if isnumeric(refSpec)
    if ~isscalar(refSpec) || refSpec < 1 || refSpec > nSessions || refSpec ~= round(refSpec)
        error('referenceSession numeric value must be an integer between 1 and %d.', nSessions);
    end
    refIdx = refSpec;
    return
end

refStr = char(string(refSpec));

% first try exact date match within datpath
refIdx = [];
for i = 1:nSessions
    datpath = char(string(out.S{i}.datpath));
    if contains(datpath, refStr)
        refIdx = i;
        break
    end
end

if isempty(refIdx)
    error('Could not resolve referenceSession "%s" to one of the loaded sessions.', refStr);
end

end