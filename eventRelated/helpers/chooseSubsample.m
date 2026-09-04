function sel = chooseSubsample(nAvail, nTake, seed)
if seed ~= 0
    rng(seed);
    sel = sort(randperm(nAvail, nTake));
else
    sel = unique(round(linspace(1, nAvail, nTake)));
    while numel(sel) < nTake
        cand = setdiff(1:nAvail, sel);
        sel = sort([sel, cand(1)]);
    end
end
end