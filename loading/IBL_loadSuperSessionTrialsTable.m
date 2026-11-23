function trialsT = IBL_loadSuperSessionTrialsTable(datpaths)
% Concatenate trials from multiple sessions, keeping session identity.

if ischar(datpaths) || isstring(datpaths)
    datpaths = cellstr(datpaths);
end

tables = cell(numel(datpaths),1);

for i = 1:numel(datpaths)

    T = IBL_loadTrialsTable(datpaths{i},'sync','timeline');

    %T = T(1:400,:); %HACK to select first 400 trials

    T.sessionID = repmat(i,height(T),1);

    tables{i} = T;
end

trialsT = vertcat(tables{:});

end