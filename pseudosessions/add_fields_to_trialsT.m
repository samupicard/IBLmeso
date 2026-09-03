function trialsT_new = add_fields_to_trialsT(trialsT)  

%adds fields to incomplete trials table (e.g. pseudoSessions)
%
% Samuel Picard

fieldnms = fieldnames(trialsT(1));
trialsT_new = trialsT;
if ~any(strcmp(fieldnms,'contrastDiff'))
    for i = 1:size(trialsT,2)
        for ii = 1:size(fieldnms,1)
            if size(trialsT(i).(fieldnms{ii}),1)==1
                trialsT_new(i).(fieldnms{ii}) = trialsT(i).(fieldnms{ii})';
            end
        end
        trialsT_new(i).contrastDiff = getContrastDiff(trialsT(i).contrastLeft,trialsT(i).contrastRight);
    end
end