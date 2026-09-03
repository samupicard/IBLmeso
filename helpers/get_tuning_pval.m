function p = get_tuning_pval(R)

%get_tuning_pval returns a p-value for a statistical test that tests if
%responses to frequency/level combinations are systematically different for
%some frequency/level combinations v. others
%
% Samuel Picard (sept 2025)

% R: nFrequencies x nLevels x nMeasurements
[nF, nL, nM] = size(R);

% Long format
[Fidx, Lidx, ~] = ndgrid(1:nF, 1:nL, 1:nM);
y = R(:);

% Make ONE grouping vector (same length as y)
% Option A: numeric group ID per (F,L)
groupID = sub2ind([nF nL], Fidx(:), Lidx(:));
Group = categorical(groupID);

% (Optional) remove NaNs
good = ~isnan(y);
y = y(good);
Group = Group(good);

% Single overall p-value: are any F×L cells different?

%parametric
%p = anova1(y, Group, 'off');  % suppress plot

%non-parametric
p = kruskalwallis(y, Group, 'off');
