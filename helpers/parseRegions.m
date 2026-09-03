function [regions, paths] = parseRegions(folder)

    files = dir(fullfile(folder, '*.mat'));

    regions = {};
    paths   = {};

    for k = 1:numel(files)
        name = files(k).name;

        tok = regexp(name, '_([^_]+)\.mat$', 'tokens');
        if isempty(tok)
            continue
        end

        regions{end+1} = tok{1}{1}; %#ok<AGROW>
        paths{end+1}   = fullfile(folder, name); %#ok<AGROW>
    end
end