function tf = hasNameValue(args, name)
% hasNameValue
%
% Return true if a name/value argument list contains the specified name.
%
% Example:
%   tf = hasNameValue(varargin, 'onlyChronic');

tf = false;

if isempty(args)
    return
end

% only inspect odd positions: {'name1', val1, 'name2', val2, ...}
for i = 1:2:numel(args)
    if ischar(args{i}) || isstring(args{i})
        if strcmpi(string(args{i}), string(name))
            tf = true;
            return
        end
    end
end
end