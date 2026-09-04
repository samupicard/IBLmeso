function [x, y, commonUID] = getPairVectors(Sodd, Seven, s, summaryMetric, eligibleUID)

uid1 = string(Sodd.chronicUID(:));
uid2 = string(Seven.chronicUID(:));

uid1 = uid1(strlength(strtrim(uid1)) > 0);
uid2 = uid2(strlength(strtrim(uid2)) > 0);

commonUID = intersect(uid1, uid2, 'stable');

if nargin >= 5 && ~isempty(eligibleUID)
    commonUID = intersect(commonUID, string(eligibleUID(:)), 'stable');
end

if isempty(commonUID)
    x = [];
    y = [];
    return
end

[row1, row2] = matchRowsByUID(Sodd.chronicUID, Seven.chronicUID, commonUID);

switch summaryMetric
    case 'diff'
        pethType = Sodd.pethTypes{s};
        [evnt, subtype] = splitPethType(pethType);

        twin1 = resolveTwinEv(Sodd.options.twin_ev, pethType, evnt, subtype);
        twin2 = resolveTwinEv(Seven.options.twin_ev, pethType, evnt, subtype);

        T1 = Sodd.T_byType{s};
        T2 = Seven.T_byType{s};

        tMask1 = T1 >= twin1(1) & T1 <= twin1(2);
        tMask2 = T2 >= twin2(1) & T2 <= twin2(2);

        x = mean(Sodd.Diff_odd{s}(row1,tMask1), 2, 'omitmissing');
        y = mean(Seven.Diff_even{s}(row2,tMask2), 2, 'omitmissing');

    case 'mi'
        x = Sodd.MI_odd{s}(row1);
        y = Seven.MI_even{s}(row2);
end

x = x(:);
y = y(:);
end