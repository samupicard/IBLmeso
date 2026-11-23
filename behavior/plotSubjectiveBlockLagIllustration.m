function plotSubjectiveBlockLagIllustration()
% plotSubjectiveBlockLagIllustration
%
% Clean paper-style illustration of objective vs subjective block identity
% for one example session. Only biased blocks (0.8 / 0.2) are shown.
%
% Positive k means the subjective prior switches later than the objective block.
%
% Example session:
%   0.8 block for 35 trials
%   0.2 block for 32 trials
%   0.8 block for 38 trials
%
% Rows shown:
%   Objective
%   k = 0
%   k = 8
%   k = 16

%% Example objective block structure
blockVals = [0.8 0.2 0.8];
blockLens = [35 32 38];
ks = [0 8 16];

% Build objective sequence
Bobj = [];
for i = 1:numel(blockVals)
    Bobj = [Bobj repmat(blockVals(i),1,blockLens(i))]; %#ok<AGROW>
end

nTrials = numel(Bobj);
x = 1:nTrials;

% Objective transition positions (between trials)
objTransitions = cumsum(blockLens(1:end-1)) + 0.5;

%% Colors
col08 = [0.72 0.42 0.05];   % brown/orange
col02 = [0.42 0.39 0.60];   % muted purple
guideCol = 0.75*[1 1 1];

%% Figure
figure('Color','w','Position',[200 200 900 500]);

% Manual axis layout for nicer spacing
nRows = 1 + numel(ks);

left = 0.12;
width = 0.83;
rowH = 0.14;

topMargin = 0.08;
bottomMargin = 0.10;

bigGap = 0.10;    % extra gap between Objective and first subjective row
smallGap = 0.04;  % regular gap between subjective rows

yTop = 1 - topMargin - rowH;

axPos = zeros(nRows,4);
axPos(1,:) = [left, yTop, width, rowH];

for iRow = 2:nRows
    if iRow == 2
        y = axPos(iRow-1,2) - bigGap - rowH;
    else
        y = axPos(iRow-1,2) - smallGap - rowH;
    end
    axPos(iRow,:) = [left, y, width, rowH];
end

%% Objective row
ax(1) = axes('Position',axPos(1,:)); %#ok<AGROW>
plotBlockRow_clean(ax(1), x, Bobj, objTransitions, col08, col02, guideCol);
%ylabel(ax(1),'Objective','FontSize',11);
%ax(1).YLabel.Rotation = 0;
%ax(1).YLabel.HorizontalAlignment = 'right';
%ax(1).YLabel.VerticalAlignment = 'middle';

%% Subjective rows
for iK = 1:numel(ks)
    k = ks(iK);

    % Shift by k trials with boundary padding
    idx = max(1, min(nTrials, (1:nTrials) - k));
    Bsubj = Bobj(idx);

    ax(iK+1) = axes('Position',axPos(iK+1,:)); %#ok<AGROW>
    plotBlockRow_clean(ax(iK+1), x, Bsubj, objTransitions, col08, col02, guideCol);
    ylabel(ax(iK+1),sprintf('k = %d',k),'FontSize',11);
    ax(iK+1).YLabel.Rotation = 0;
    ax(iK+1).YLabel.HorizontalAlignment = 'right';
    ax(iK+1).YLabel.VerticalAlignment = 'middle';
end

%% Shared axis cosmetics
for i = 1:nRows
    if i < nRows
        set(ax(i),'XTickLabel',[]);
    else
        xlabel(ax(i),'Trial','FontSize',11);
    end
end

linkaxes(ax,'x');

end


function plotBlockRow_clean(ax, x, B, objTransitions, col08, col02, guideCol)
% Plot one row of colored blocks in a clean style

axes(ax); %#ok<LAXES>
hold(ax,'on');

% Find contiguous segments in this row
edges = [1, find(diff(B)~=0)+1, numel(B)+1];

for iSeg = 1:numel(edges)-1
    i1 = edges(iSeg);
    i2 = edges(iSeg+1)-1;

    if abs(B(i1)-0.8) < 1e-9
        thisColor = col08;
    else
        thisColor = col02;
    end

    patch(ax, ...
        [x(i1) x(i2) x(i2) x(i1)], ...
        [0 0 1 1], ...
        thisColor, ...
        'EdgeColor','none');

    % Optional text label inside segment
    xc = mean([x(i1) x(i2)]);
    text(ax, xc, 0.5, sprintf('%.1f',B(i1)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color','w', ...
        'FontWeight','bold', ...
        'FontSize',10);
end

% Draw objective switch-point guide lines on every row
for t = objTransitions
    xline(ax, t, ':', 'Color', guideCol, 'LineWidth', 1.2);
end

% Cosmetics
xlim(ax,[1 x(end)]);
ylim(ax,[0 1]);

ax.YTick = [];
ax.Box = 'off';
ax.TickDir = 'out';
ax.FontSize = 10;
ax.LineWidth = 1;
ax.XColor = [0.2 0.2 0.2];
ax.YColor = [0.2 0.2 0.2];

end