function foldId = makeContiguousFolds(T, nFolds)
%MAKECONTIGUOUSFOLDS Divide observations into contiguous temporal folds.

if nFolds < 3
    error("makeContiguousFolds:TooFewFolds", ...
        "At least three folds are required.");
end

if nFolds > T
    error("makeContiguousFolds:TooManyFolds", ...
        "nFolds cannot exceed the number of observations.");
end

foldEdges = round(linspace(0, T, nFolds + 1));
foldId = zeros(T,1);

for foldIndex = 1:nFolds
    firstIndex = foldEdges(foldIndex) + 1;
    lastIndex = foldEdges(foldIndex + 1);

    foldId(firstIndex:lastIndex) = foldIndex;
end

end