function [regionIdsExpanded, regsLabel] = expand_region_ids(opt)
% Returns:
%   regionIdsExpanded: vector of Allen IDs to keep (empty means "no filter" or "no matches")
%   regsLabel: label string for plots/logs

regionIdsExpanded = [];

% normalize opt.region to cellstr
if isempty(opt) || ~isfield(opt,'region') || isempty(opt.region)
    regsLabel = "All";
    return
end

if isstring(opt.region), opt.region = cellstr(opt.region); end
if ischar(opt.region),   opt.region = {opt.region}; end

% get structure tree
st = [];
if isfield(opt,'st') && ~isempty(opt.st)
    st = opt.st;
elseif isfield(opt,'stPath') && strlength(string(opt.stPath)) > 0
    st = loadStructureTree(opt.stPath);
else
    error('Region filter requested but no structure tree provided. Set opt.st (table) or opt.stPath.');
end

% normalize fields
acr_all = string(st.acronym);
rid_all = double(st.id);

if ismember('depth', st.Properties.VariableNames)
    depth_all = double(st.depth);
else
    error('Structure tree table must contain a "depth" column.');
end

% 1) expand tokens (e.g. HVAs -> {'VISa','VISl',...}) and keep others as-is
regList = {};  % cellstr
for r = 1:numel(opt.region)
    tok = opt.region{r};

    acrList = regionTokenMap(tok);
    if ~isempty(acrList)
        regList = [regList, cellstr(string(acrList(:)))']; %#ok<AGROW>
    else
        regList = [regList, {tok}]; %#ok<AGROW>
    end
end

% 2) build union of matching IDs
for r = 1:numel(regList)
    reg = char(regList{r});

    % wildcard: 'VIS*'
    if endsWith(reg, '*')
        prefix = reg(1:end-1);
        match = startsWith(acr_all, string(prefix));
        regionIdsExpanded = [regionIdsExpanded; rid_all(match)]; %#ok<AGROW>
        continue
    end

    % exact acronym: 'VISa'
    exact = (acr_all == string(reg));
    if ~any(exact)
        continue
    end

    regionIdsExpanded = [regionIdsExpanded; rid_all(exact)]; %#ok<AGROW>

    % depth-gated prefix expansion if depth >= 7
    d = depth_all(find(exact, 1, 'first'));
    if d >= 7
        match = startsWith(acr_all, string(reg));
        regionIdsExpanded = [regionIdsExpanded; rid_all(match)]; %#ok<AGROW>
    end
end

regionIdsExpanded = unique(regionIdsExpanded);

% label
if numel(opt.region) == 1
    regsLabel = string(opt.region{1});
else
    regsLabel = string(strjoin(opt.region, ', '));
end
end
