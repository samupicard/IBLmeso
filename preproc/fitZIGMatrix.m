function zig = fitZIGMatrix(X)

nROIs = size(X,2);

zig.q   = nan(nROIs,1);
zig.loc = nan(nROIs,1);
zig.k   = nan(nROIs,1);
zig.a   = nan(nROIs,1);
zig.nPos = nan(nROIs,1);

for iROI = 1:nROIs
    z = fitZeroInflatedGamma(X(:,iROI));

    zig.q(iROI)    = z.q;
    zig.loc(iROI)  = z.loc;
    zig.k(iROI)    = z.k;
    zig.a(iROI)    = z.a;
    zig.nPos(iROI) = z.nPos;
end

end