function [row1, row2] = matchRowsByUID(uidVec1, uidVec2, uidKeep)

uidVec1 = string(uidVec1(:));
uidVec2 = string(uidVec2(:));
uidKeep = string(uidKeep(:));

row1 = nan(numel(uidKeep),1);
row2 = nan(numel(uidKeep),1);

for i = 1:numel(uidKeep)
    row1(i) = find(uidVec1 == uidKeep(i), 1, 'first');
    row2(i) = find(uidVec2 == uidKeep(i), 1, 'first');
end

good = isfinite(row1) & isfinite(row2);
row1 = row1(good);
row2 = row2(good);
end