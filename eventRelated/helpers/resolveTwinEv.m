function twin = resolveTwinEv(twinOpt, pethType, evnt, subtype)
if isempty(twinOpt)
    [defaultTwinEv, ~, ~] = get_twins(evnt, subtype);
    twin = defaultTwinEv;
    return
end

if isnumeric(twinOpt)
    twin = twinOpt;
    return
end

if isstruct(twinOpt)
    fn = matlab.lang.makeValidName(pethType);
    if isfield(twinOpt, fn)
        twin = twinOpt.(fn);
    else
        [defaultTwinEv, ~, ~] = get_twins(evnt, subtype);
        twin = defaultTwinEv;
    end
    return
end

twin = [];
end