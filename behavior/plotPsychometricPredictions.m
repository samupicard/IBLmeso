function plotPsychometricPredictions(trialsT, results, varargin)
%PLOTPSYCHOMETRICPREDICTIONS Plot observed and predicted choices across trials.
%
% plotPsychometricPredictions(trialsT,results)
% plotPsychometricPredictions(trialsT,results,'cv',false)
% plotPsychometricPredictions(trialsT,results,'cv',true,'folds',1:3)
% plotPsychometricPredictions(trialsT,results,'cv',false,'maxTrials',500)
%
% NAME-VALUE OPTIONS
%   'cv'        true (default): out-of-fold predictions and fold parameters;
%               false: predictions reconstructed from results.full.
%   'folds'     'all' (default), or numeric CV fold IDs; only for cv=true.
%   'maxTrials' [] (default), or maximum contiguous original trials shown.
%
% In CV mode, each held-out trial uses the model and bias parameters of its
% respective training fold. In full-fit mode, models use the coefficients
% and parameters in results.full, including results.full.drift.biasPredict.
% Only results.validTrials are scored/marked in either mode.
%
% All five model traces show the total fitted bias in log-odds units.
% Incorrectly classified trials (threshold 0.5) are marked by crosses,
% coloured according to the ACTUAL choice. An overall error rate is shown.
%
% Assumes the model/history conventions used by comparePsychometricModelsCV.

p = inputParser;
p.addParameter('cv',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('folds','all',@(x)(ischar(x)||isstring(x)||isnumeric(x)));
p.addParameter('maxTrials',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1));
p.parse(varargin{:});
opt = p.Results;

n = height(trialsT);
assert(numel(results.validTrials)==n,'results.validTrials must have one entry per trial.');
valid = logical(results.validTrials(:));

% Choice convention: right = +1 after IBL sign inversion.
choice = -double(trialsT.choice(:));
actualRightAll = choice==1;

if ismember('contrastDiff',trialsT.Properties.VariableNames)
    X = double(trialsT.contrastDiff(:));
else
    cr = double(trialsT.contrastRight(:));
    cl = double(trialsT.contrastLeft(:));
    cr(isnan(cr)) = 0;
    cl(isnan(cl)) = 0;
    X = cr-cl;
end
pLeft = double(trialsT.probabilityLeft(:));

% Determine boundaries independently for each session in a superSession.
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

% Convention: +1 right-favouring block, -1 left-favouring, 0 unbiased.
B = zeros(n,1);
B(abs(pLeft-.2)<1e-6) =  1;
B(abs(pLeft-.8)<1e-6) = -1;

% Regressor/prediction arrays, populated at appropriate eligible indices.
% Columns: baseline, block, stimulus side, contrast-weighted, drift.
pred = nan(n,5);
bias = nan(n,5);
fX = nan(n,1);

if opt.cv
    foldId = results.foldId(:);
    assert(numel(foldId)==n,'results.foldId must have one entry per trial.');
    if (ischar(opt.folds)||isstring(opt.folds)) && strcmpi(string(opt.folds),'all')
        folds = unique(foldId(valid & foldId>0))';
    elseif isnumeric(opt.folds)
        folds = unique(opt.folds(:)');
    else
        error('folds must be numeric or ''all''.');
    end
    selected = valid & ismember(foldId,folds);
    assert(any(selected),'No eligible held-out trials for selected folds.');
    assert(isfield(results,'foldBeta') && isfield(results,'cvPrediction') && ...
        isfield(results,'cvBiasPredict'), ...
        'CV mode requires foldBeta, cvPrediction and cvBiasPredict in results.');
    assert(size(results.cvPrediction,2)==5, ...
        'Expected five columns in results.cvPrediction.');
    pred(selected,:) = results.cvPrediction(selected,:);

    for f = folds
        idx = selected & foldId==f;
        if ~any(idx), continue; end
        q = results.folds.q(f);
        k = results.folds.bestK(f);
        tauSide = results.folds.bestTauSide(f);
        tauContrast = results.folds.bestTauContrast(f);
        transformed = sign(X).*abs(X).^q;
        fX(idx) = transformed(idx);
        Bshift = shiftBlock(B,k,starts,ends);
        Hside = stimulusHistory(X,q,tauSide,starts,ends,false);
        Hcontrast = stimulusHistory(X,q,tauContrast,starts,ends,true);
        base = results.foldBeta{f,1};
        block = results.foldBeta{f,2};
        side = results.foldBeta{f,3};
        contrast = results.foldBeta{f,4};
        drift = results.cvBiasPredict(:);
        bias(idx,1) = base(1);
        bias(idx,2) = block(1)+block(3)*Bshift(idx);
        bias(idx,3) = side(1)+side(3)*Hside(idx);
        bias(idx,4) = contrast(1)+contrast(3)*Hcontrast(idx);
        bias(idx,5) = base(1)+drift(idx);
    end
    modeLabel = sprintf('Held-out predictions | folds %s',mat2str(folds));
else
    selected = valid;
    assert(any(selected),'No eligible trials in results.validTrials.');
    full = results.full;
    % These are the fields saved by comparePsychometricModelsCV.
    q = full.q;
    k = full.block.bestK;
    tauSide = full.side.bestTau;
    tauContrast = full.contrast.bestTau;
    base = full.baseline.beta(:);
    block = full.block.beta(:);
    side = full.side.beta(:);
    contrast = full.contrast.beta(:);
    assert(isfield(full.drift,'biasPredict'), ...
        'Full drifting model must contain full.drift.biasPredict.');
    drift = double(full.drift.biasPredict(:));
    assert(numel(drift)==n, ...
        'results.full.drift.biasPredict must have one element per original trial.');
    transformed = sign(X).*abs(X).^q;
    fX(selected) = transformed(selected);
    Bshift = shiftBlock(B,k,starts,ends);
    Hside = stimulusHistory(X,q,tauSide,starts,ends,false);
    Hcontrast = stimulusHistory(X,q,tauContrast,starts,ends,true);

    % Total fitted bias, excluding contrast-dependent stimulus term.
    bias(selected,1) = base(1);
    bias(selected,2) = block(1)+block(3)*Bshift(selected);
    bias(selected,3) = side(1)+side(3)*Hside(selected);
    bias(selected,4) = contrast(1)+contrast(3)*Hcontrast(selected);
    bias(selected,5) = base(1)+drift(selected);

    % Reconstruct logistic probabilities from the complete-data fits.
    % Drift inherits stimulus coefficients from the baseline model.
    pred(selected,1) = logistic(base(1)+base(2)*transformed(selected));
    pred(selected,2) = logistic(block(1)+block(2)*transformed(selected)+block(3)*Bshift(selected));
    pred(selected,3) = logistic(side(1)+side(2)*transformed(selected)+side(3)*Hside(selected));
    pred(selected,4) = logistic(contrast(1)+contrast(2)*transformed(selected)+contrast(3)*Hcontrast(selected));
    % The drift model uses a latent random-walk bias; this is its
    % one-step-ahead *state estimate* plugged into a logistic link.
    % If full.drift stores integrated predictive probabilities, prefer
    % those for exact comparison to the original quadrature model.
    pred(selected,5) = logistic(base(1)+base(2)*transformed(selected)+drift(selected));
    modeLabel = 'Full-data model predictions (in-sample)';
end

% Restrict to a contiguous original trial interval; retain gaps in markers.
selIdx = find(selected);
iStart = selIdx(1);
iEnd = selIdx(end);
if ~isempty(opt.maxTrials)
    iEnd = min(iEnd,iStart+floor(opt.maxTrials)-1);
end
plotIdx = (iStart:iEnd)';
xChoice = plotIdx(selected(plotIdx));
yStim = fX(xChoice);
actualRight = actualRightAll(xChoice);

% Shared axes.
contrastMax = max(abs(fX(xChoice)),[],'omitnan');
if isempty(contrastMax)||~isfinite(contrastMax)||contrastMax==0
    contrastMax = 1;
end
contrastYLim = 1.1*[-contrastMax contrastMax];
allBias = bias(plotIdx,:);
biasLim = max(abs(allBias(:)),[],'omitnan');
if isempty(biasLim)||~isfinite(biasLim), biasLim=0; end
biasLim = 1.1*max(.5,biasLim);
biasYLim = [-biasLim biasLim];

colRight = [.2 .4 .8];
colLeft = [.8 .25 .2];
colBias = [.12 .48 .35];
col02 = [.95 .78 .49];
col08 = [.68 .62 .82];
actualColors = repmat(colLeft,numel(xChoice),1);
actualColors(actualRight,:) = repmat(colRight,sum(actualRight),1);

modelNames = {'Actual choices','Baseline','Block', ...
    'Stimulus side','Contrast weighted','Drifting bias'};

figure('Color','w','Position',[150 80 1250 850]);
tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');
ax = gobjects(6,1);

for iModel = 1:6
    ax(iModel) = nexttile(tl);
    hold(ax(iModel),'on');
    yyaxis(ax(iModel),'left');
    ylim(ax(iModel),contrastYLim);
    ax(iModel).YAxis(1).Color = [.2 .2 .2];

    % Subtle objective-block backgrounds on contrast (left) axis.
    blockChanges = [1;find(diff(pLeft)~=0)+1;n+1];
    for j = 1:numel(blockChanges)-1
        a = max(blockChanges(j)-0.5,iStart-0.5);
        b = min(blockChanges(j+1)-0.5,iEnd+0.5);
        if b<=a,continue;end
        prob = pLeft(blockChanges(j));
        if abs(prob-.2)<1e-6
            col = col02;
        elseif abs(prob-.8)<1e-6
            col = col08;
        else
            continue
        end
        patch(ax(iModel),[a b b a], ...
            [contrastYLim(1) contrastYLim(1) contrastYLim(2) contrastYLim(2)], ...
            col,'EdgeColor','none','FaceAlpha',.20, ...
            'HandleVisibility','off');
    end

    if iModel>1
        yyaxis(ax(iModel),'right');
        % Separate each fold's bias curve: don't connect independently
        % fitted models, or bridge non-selected trials.
        if opt.cv
            for f = folds
                ii = plotIdx(selected(plotIdx) & foldId(plotIdx)==f & ...
                    isfinite(bias(plotIdx,iModel-1)));
                plot(ax(iModel),ii,bias(ii,iModel-1),'-', ...
                    'Color',colBias,'LineWidth',1.2);
            end
        else
            ii = plotIdx(selected(plotIdx) & isfinite(bias(plotIdx,iModel-1)));
            plot(ax(iModel),ii,bias(ii,iModel-1),'-', ...
                'Color',colBias,'LineWidth',1.2);
        end
        ylim(ax(iModel),biasYLim);
        ax(iModel).YAxis(2).Color = colBias;
        ylabel(ax(iModel),'Bias (log-odds)');
        yyaxis(ax(iModel),'left');
    end

    if iModel==1
        scatter(ax(iModel),xChoice,yStim,12,actualColors);
    else
        pPred = pred(xChoice,iModel-1);
        hasPred = isfinite(pPred);
        predictedRight = pPred>=.5;
        correct = predictedRight==actualRight;
        correct = correct & hasPred;
        wrong = ~correct & hasPred;
        brierScore = mean((pPred(hasPred)-double(actualRight(hasPred))).^2);
        % predictedColors = repmat(colLeft,numel(xChoice),1);
        % predictedColors(predictedRight,:) = repmat(colRight,sum(predictedRight),1);
        % scatter(ax(iModel),xChoice(correct),yStim(correct), ...
        %     12,predictedColors(correct,:),'filled');
        scatter(ax(iModel),xChoice(wrong),yStim(wrong), ...
            20,actualColors(wrong,:),'x','LineWidth',1.5);
        if any(hasPred)
            pWrong = mean(wrong(hasPred));
            text(ax(iModel),.98,.95,sprintf('Error: %.1f%%',100*pWrong), ...
                'Units','normalized','HorizontalAlignment','right', ...
                'VerticalAlignment','top','FontSize',9, ...
                'Color',[.25 .25 .25]);
            text(ax(iModel),.02,.95,sprintf('Brier: %.3f',brierScore), ...
                'Units','normalized', ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','top', ...
                'FontSize',9, ...
                'Color',[.25 .25 .25]);
        end
    end

    % Session boundaries, shared across both plotting modes.
    for j = 2:numel(starts)
        if starts(j)>iStart && starts(j)<=iEnd
            xline(ax(iModel),starts(j)-.5,':', ...
                'Color',[.45 .45 .45],'HandleVisibility','off');
        end
    end
    if opt.cv
        for j = 2:numel(plotIdx)
            t0 = plotIdx(j-1);t1 = plotIdx(j);
            if foldId(t0)>0 && foldId(t1)>0 && foldId(t0)~=foldId(t1)
                xline(ax(iModel),t1-.5,'--', ...
                    'Color',[.25 .55 .25], ...
                    'HandleVisibility','off');
            end
        end
    end

    ylim(ax(iModel),contrastYLim);
    yticks(ax(iModel),[-contrastMax 0 contrastMax]);
    yticklabels(ax(iModel),{'Left','0','Right'});
    ylabel(ax(iModel),modelNames{iModel});
    xlim(ax(iModel),[iStart-.5 iEnd+.5]);
    ax(iModel).TickDir = 'out';
    ax(iModel).FontSize = 9;
    box(ax(iModel),'off');
    if iModel<6,ax(iModel).XTickLabel=[];end
end
xlabel(ax(end),'Trial');
linkaxes(ax,'x');
title(tl,modeLabel);
end

function Bshift = shiftBlock(B,k,starts,ends)
Bshift = zeros(size(B));
for s = 1:numel(starts)
    ids = (starts(s):ends(s))';
    idx = max(1,min(numel(ids),(1:numel(ids))'-k));
    Bshift(ids) = B(ids(idx));
end
end

function H = stimulusHistory(X,q,tau,starts,ends,weighted)
% Exponentially weighted mean over preceding nonzero-contrast stimuli.
% Match the helper in comparePsychometricModelsCV if it differs.
H = zeros(size(X));
decay = exp(-1/tau);
for s = 1:numel(starts)
    numerator=0;denominator=0;
    for t = starts(s):ends(s)
        if denominator>0,H(t)=numerator/denominator;end
        if ~isfinite(X(t))||X(t)==0,continue;end
        if weighted,w=abs(X(t))^q;else,w=1;end
        numerator=decay*numerator+w*sign(X(t));
        denominator=decay*denominator+w;
    end
end
end

function p = logistic(z)
p = 1./(1+exp(-z));
end
