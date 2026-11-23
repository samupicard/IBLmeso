function plotPsychometricCVPredictions(trialsT,results,fold,maxTrials)
%PLOTPSYCHOMETRICCVPREDICTIONS Visualize held-out model predictions.
%
% plotPsychometricCVPredictions(trialsT,results,fold,maxTrials)
%
% Rows:
%   1. Actual choices
%   2. Baseline
%   3. Discrete subjective block
%   4. Stimulus-side history
%   5. Contrast-weighted history
%   6. Drifting bias
%
% Markers:
%   Y-position : signed stimulus contrast, transformed by exponent q
%   Blue       : right choice
%   Red        : left choice
%   Filled     : correct prediction
%   Cross      : incorrect prediction
%
% Marker colours reflect actual choices in the first row, and
% predicted choices in all subsequent rows.
%
% Lines:
%   Green      : total bias contribution (log-odds)
%   Shading    : objective probabilityLeft blocks
%
% All contrast axes share the same limits.
% All bias axes share the same limits, independently scaled.

if nargin<3 || isempty(fold), fold = 1; end
if nargin<4, maxTrials = []; end

%% Extract trial variables

n = height(trialsT);

if ismember('contrastDiff',trialsT.Properties.VariableNames)
    X = double(trialsT.contrastDiff(:));
else
    cr = double(trialsT.contrastRight(:));
    cl = double(trialsT.contrastLeft(:));
    cr(isnan(cr)) = 0;
    cl(isnan(cl)) = 0;
    X = cr-cl;
end

choice = -double(trialsT.choice(:));
y = double(choice==1);
pLeft = double(trialsT.probabilityLeft(:));

%% Select folds and display range

if ischar(fold) || isstring(fold)
    if strcmpi(fold,'all')
        fold = 1:max(results.foldId);
    else
        error('Unknown fold selection.');
    end
end

fold = unique(fold(:)');

valid = results.validTrials(:);
test = valid & ismember(results.foldId(:),fold);

testIdx = find(test);
assert(~isempty(testIdx),'No test trials for selected folds.');

% Contiguous segment of the original trial sequence
iStart = testIdx(1);
iEnd = testIdx(end);

if ~isempty(maxTrials)
    iEnd = min(iEnd,iStart+maxTrials-1);
end

plotIdx = (iStart:iEnd)';
xPlot = plotIdx;

% Held-out trials in the displayed range
xChoice = plotIdx(test(plotIdx));
yChoice = y(xChoice);

% Fold identity for each plotted choice
choiceFold = results.foldId(xChoice);

%% Identify session boundaries

if ismember('sessionID',trialsT.Properties.VariableNames)
    sessionCode = findgroups(trialsT.sessionID);
elseif ismember('sessionNr',trialsT.Properties.VariableNames)
    sessionCode = findgroups(trialsT.sessionNr);
else
    sessionCode = ones(n,1);
end

sessionCode = sessionCode(:);

starts = [1;find(diff(sessionCode)~=0)+1];
ends = [starts(2:end)-1;n];

%% Reconstruct fold-specific bias contributions

% Matrix of total bias contributions:
% columns = baseline, block, side, contrast, drift
biasMatrix = nan(n,5);

% Transformed contrast used for marker y-position
fX = nan(n,1);

% Loop over requested folds
for iFold = 1:numel(fold)

    f = fold(iFold);

    % Held-out trials belonging to this fold
    testF = valid & results.foldId(:)==f;

    %% Fold-specific parameters

    q = results.folds.q(f);
    k = results.folds.bestK(f);

    tauSide = results.folds.bestTauSide(f);
    tauContrast = results.folds.bestTauContrast(f);

    fXf = sign(X).*abs(X).^q;
    fX(testF) = fXf(testF);

    %% Subjective block regressor

    B = zeros(n,1);
    B(abs(pLeft-0.2)<1e-6) = 1;
    B(abs(pLeft-0.8)<1e-6) = -1;

    Bshift = zeros(n,1);

    for s = 1:numel(starts)
        ids = (starts(s):ends(s))';
        idx = max(1,min(numel(ids),(1:numel(ids))'-k));
        Bshift(ids) = B(ids(idx));
    end

    %% Stimulus-history regressors

    Hside = makeHistory(X,tauSide,q,starts,ends,false);
    Hcontrast = makeHistory(X,tauContrast,q,starts,ends,true);

    %% Fitted coefficients

    betaBase = results.foldBeta{f,1};
    betaBlock = results.foldBeta{f,2};
    betaSide = results.foldBeta{f,3};
    betaContrast = results.foldBeta{f,4};

    %% Total bias contributions on held-out trials

    biasMatrix(testF,1) = betaBase(1);

    biasMatrix(testF,2) = ...
        betaBlock(1) + betaBlock(3)*Bshift(testF);

    biasMatrix(testF,3) = ...
        betaSide(1) + betaSide(3)*Hside(testF);

    biasMatrix(testF,4) = ...
        betaContrast(1) + betaContrast(3)*Hcontrast(testF);

    biasMatrix(testF,5) = ...
        betaBase(1) + results.cvBiasPredict(testF);

end

%% Bias traces for the six plotting rows

biasTraces = { ...
    [], ...
    biasMatrix(:,1), ...
    biasMatrix(:,2), ...
    biasMatrix(:,3), ...
    biasMatrix(:,4), ...
    biasMatrix(:,5)};


%% Model labels

modelNames = { ...
    'Actual choices', ...
    'Baseline', ...
    'Block', ...
    'Stimulus side', ...
    'Contrast weighted', ...
    'Drifting bias'};

nRows = numel(modelNames);

%% Colours

colRight = [0.20 0.40 0.80];  % blue
colLeft  = [0.80 0.25 0.20];  % red

colBias = [0.12 0.48 0.35];   % dark muted green

col02 = [0.95 0.78 0.49];     % orange block
col08 = [0.68 0.62 0.82];     % purple block

%% Common axis limits

% Contrast axis: fixed by the transformed stimulus range
contrastMax = max(abs(fX(xChoice)),[],'omitnan');

if isempty(contrastMax) || contrastMax==0
    contrastMax = 1;
end

contrastLim = contrastMax*1.10;
contrastYLim = [-contrastLim contrastLim];

% Bias axis: fixed across ALL model rows
allBias = [];

for iModel = 2:nRows
    b = biasTraces{iModel}(plotIdx);
    allBias = [allBias; b(isfinite(b))]; %#ok<AGROW>
end

if isempty(allBias)
    biasLim = 1;
else
    biasLim = max(0.5,1.10*max(abs(allBias)));
end

biasYLim = [-biasLim biasLim];

% The two axes have different units but share the same
% physical vertical extent and zero position.

%% Prepare choice marker positions and colours

yStim = fX(xChoice);

% Actual choices
actualRight = logical(yChoice);

actualColors = repmat(colLeft,numel(xChoice),1);
actualColors(actualRight,:) = ...
    repmat(colRight,sum(actualRight),1);

%% Figure

figure('Color','w','Position',[100 100 2000 850]);

tl = tiledlayout(nRows,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax = gobjects(nRows,1);

for iModel = 1:nRows

    ax(iModel) = nexttile;
    hold(ax(iModel),'on');

    %% Left axis: contrast and block backgrounds
    if iModel>1
        yyaxis(ax(iModel),'left');
    end

    ylim(ax(iModel),contrastYLim);
    yticks(ax(iModel),[-contrastMax 0 contrastMax]);

    if iModel==1
        yticklabels(ax(iModel),{'Left','0','Right'});
    else
        yticklabels(ax(iModel),[]);
    end

    ax(iModel).YAxis(1).Color = [0.25 0.25 0.25];

    % Zero-contrast reference
    yline(ax(iModel),0,':', ...
        'Color',[0.6 0.6 0.6], ...
        'HandleVisibility','off');

    % Objective block shading
    thisP = pLeft(plotIdx);

    edges = [1;find(diff(thisP)~=0)+1;numel(plotIdx)+1];

    for j = 1:numel(edges)-1

        i1 = edges(j);
        i2 = edges(j+1)-1;

        if abs(thisP(i1)-0.2)<1e-6
            col = col02;
        elseif abs(thisP(i1)-0.8)<1e-6
            col = col08;
        else
            continue
        end

        patch(ax(iModel), ...
            [plotIdx(i1)-0.5 plotIdx(i2)+0.5 ...
            plotIdx(i2)+0.5 plotIdx(i1)-0.5], ...
            [contrastYLim(1) contrastYLim(1) ...
            contrastYLim(2) contrastYLim(2)], ...
            col, ...
            'FaceAlpha',0.18, ...
            'EdgeColor','none', ...
            'HandleVisibility','off');
    end

    %% Right axis: total model bias

    if iModel>1

        yyaxis(ax(iModel),'right');

        thisBias = biasTraces{iModel};

        for f = fold
            ii = plotIdx( ...
                results.foldId(plotIdx)==f & isfinite(thisBias(plotIdx)));

            if ~isempty(ii)
                plot(ax(iModel),ii,thisBias(ii),'-', ...
                    'Color',colBias,'LineWidth',1.5);
            end
        end

        ylim(ax(iModel),biasYLim);

        % Colour the right axis to match the bias trace
        ax(iModel).YAxis(2).Color = colBias;

        ylabel(ax(iModel),'Bias (log-odds)');
        %ax(iModel).YAxis(2).TickLabels = {};

        % Return to contrast axis for markers
        yyaxis(ax(iModel),'left');
    end

    %% Plot choices
    if iModel==1

        % Ground truth: actual choices
        scatter(ax(iModel),xChoice,yStim,12,actualColors, ...
            'LineWidth',1);%'filled','MarkerEdgeColor','none');

    else
        if false

            % Model predictions on held-out trials
            pPred = results.cvPrediction(xChoice,iModel-1);

            predictedRight = pPred>=0.5;
            correct = predictedRight==actualRight;

            % % Colour markers by predicted choice
            % colors = repmat(colLeft,numel(xChoice),1);
            % colors(predictedRight,:) = ...
            %     repmat(colRight,sum(predictedRight),1);

            % % Correct predictions: filled circles
            % scatter(ax(iModel), ...
            %     xChoice(correct),yStim(correct), ...
            %     20,colors(correct,:), ...
            %     'filled','MarkerEdgeColor','none');

            % Incorrect predictions: crosses
            scatter(ax(iModel), ...
                xChoice(~correct),yStim(~correct), ...
                25,actualColors(~correct,:), ...
                'x','LineWidth',1.6);

            % Compute fractions of incorrect predictions
            pWrongLeft  = sum(~correct & ~actualRight) / sum(~actualRight);
            pWrongRight = sum(~correct & actualRight) / sum(actualRight);
            pWrong = mean(~correct);

            % Display error rates (using the corresponding choice colours)
            text(ax(iModel),0.98,0.95, ...
                sprintf('L error: %.1f%%',100*pWrongLeft), ...
                'Units','normalized', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',9,'Color',colLeft);

            text(ax(iModel),0.98,0.78, ...
                sprintf('R error: %.1f%%',100*pWrongRight), ...
                'Units','normalized', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',9,'Color',colRight);

            text(ax(iModel),0.005,0.95, ...
                sprintf('Error: %.1f%%',100*pWrong), ...
                'Units','normalized', ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','top', ...
                'FontSize',9, ...
                'Color',[0.25 0.25 0.25]);

        end
    end

    %% Session boundaries

    for s = 2:numel(starts)

        boundary = starts(s)-0.5;

        if boundary>=plotIdx(1) && boundary<=plotIdx(end)
            xline(ax(iModel),boundary,':', ...
                'Color',[0.45 0.45 0.45], ...
                'LineWidth',0.8, ...
                'HandleVisibility','off');
        end
    end

    %% CV fold boundaries

    foldIDs = results.foldId(:);

    for i = 2:numel(plotIdx)

        t1 = plotIdx(i-1);
        t2 = plotIdx(i);

        % Only mark transitions between selected, nonzero folds
        if foldIDs(t1)>0 && foldIDs(t2)>0 && ...
                foldIDs(t1)~=foldIDs(t2)

            xline(ax(iModel),t2-0.5,'--', ...
                'Color',[0.2 0.6 0.2], ...
                'LineWidth',1, ...
                'HandleVisibility','off');
        end
    end

    %% Cosmetics

    xlim(ax(iModel),[plotIdx(1)-0.5 plotIdx(end)+0.5]);

    % Ensure contrast is the active axis
    if iModel>1
        yyaxis(ax(iModel),'left');
    end
    ylim(ax(iModel),contrastYLim);

    ylabel(ax(iModel),modelNames{iModel}, ...
        'Rotation',0, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'Color',[0.2 0.2 0.2]);

    box(ax(iModel),'off');
    set(ax(iModel),'TickDir','out','Layer','top');

    if iModel<nRows
        ax(iModel).XTickLabel = [];
    else
        xlabel(ax(iModel),'Trial');
    end

end

linkaxes(ax,'x');

title(tl,sprintf('Held-out predictions, folds %s',mat2str(fold)));

end


function H = makeHistory(X,tau,q,starts,ends,weighted)
% Exponentially weighted preceding nonzero stimulus history.
% Resets at session boundaries.

alpha = exp(-1/tau);
H = zeros(size(X));

for s = 1:numel(starts)

    numerator = 0;
    denominator = 0;

    for t = starts(s):ends(s)

        if denominator>0
            H(t) = numerator/denominator;
        end

        x = X(t);

        if isfinite(x) && x~=0

            if weighted
                w = abs(x)^q;
            else
                w = 1;
            end

            numerator = alpha*numerator + w*sign(x);
            denominator = alpha*denominator + w;

        end
    end
end

end
