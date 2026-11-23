function [fig, out] = plotBlockTransitionBias(trialsT, results, varargin)
%PLOTBLOCKTRANSITIONBIAS Plot full-fit bias around block switches.
%
% [fig,out] = plotBlockTransitionBias(trialsT,results)
% [fig,out] = plotBlockTransitionBias(...,'window',[-30 60])
% [fig,out] = plotBlockTransitionBias(...,'quantity','regressor')
%
% Columns: P(left) 0.2 -> 0.8, and 0.8 -> 0.2.
% Rows: delayed block, stimulus-side history, contrast-weighted history,
%       drifting bias.
% One line per session, averaged over that session's transitions of the
% corresponding type. Colour identifies the session, consistently in all
% eight panels. Trial zero is the first trial in the new block.
%
% By default the plotted quantity is the FULL-FIT TOTAL BIAS CONTRIBUTION
% (intercept + fitted coefficient * regressor) in log-odds, as in
% plotPsychometricPredictions(...,'cv',false). Use 'quantity','regressor'
% for raw regressors (note that the drift regressor is in log-odds units).
%
% NAME-VALUE INPUTS
%   'window'     [before after] in trials, default [-10 30]
%   'quantity'   'bias' (default) or 'regressor'
%   'minTrials'  Minimum valid values per session trace to show (default 2)
%
% OUTPUTS
%   fig          Figure handle
%   out          Struct with .relativeTrials, .sessionLabels, .colors,
%                .traces (nRelativeTrials x nSessions x 4 x 2),
%                .nTransitions (nSessions x 2), .transitionIndices.
%
% Important: Uses ALL trials for calculating histories and extracting
% transition windows, not only the trials included in model fitting.
% The drifting bias is the one-step-ahead latent bias prediction saved in
% results.full.drift.biasPredict (not the posterior using current choice).

ip = inputParser;
ip.addParameter('window',[-10 30],@(x)isnumeric(x)&&numel(x)==2&& ...
    all(isfinite(x))&&all(x==round(x))&&x(1)<=0&&x(2)>=0);
ip.addParameter('quantity','bias',@(x)any(strcmpi(string(x),["bias","regressor"])));
ip.addParameter('minTrials',2,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('sessionEdges',[],@(x)isnumeric(x)&&isvector(x)&& ...
    all(isfinite(x))&&all(x==round(x)));
ip.parse(varargin{:});
opt = ip.Results;
rel = (opt.window(1):opt.window(2))';

n = height(trialsT);
assert(n>0,'trialsT is empty.');
assert(isfield(results,'full'),'results.full is required.');
full = results.full;

% Signed contrast (right minus left).
if ismember('contrastDiff',trialsT.Properties.VariableNames)
    X = double(trialsT.contrastDiff(:));
else
    cr = double(trialsT.contrastRight(:));
    cl = double(trialsT.contrastLeft(:));
    cr(isnan(cr))=0; cl(isnan(cl))=0;
    X = cr-cl;
end
pLeft = double(trialsT.probabilityLeft(:));

% Preserve the order in which session IDs occur.
if ismember('sessionID',trialsT.Properties.VariableNames)
    id = string(trialsT.sessionID(:));
elseif ismember('sessionNr',trialsT.Properties.VariableNames)
    id = string(trialsT.sessionNr(:));
else
    id = repmat("Session 1",n,1);
end
assert(all(~ismissing(id)),'Missing session ID in trialsT.');
starts = [1;find(id(2:end)~=id(1:end-1))+1];
ends = [starts(2:end)-1;n];
labels = id(starts);
nSessions = numel(starts);
if numel(unique(labels))~=nSessions
    % Ensure distinct colours even when ID values are reused in separate runs.
    labels = labels + " (segment " + string((1:nSessions)') + ")";
end

% Objective block: +1 = P(left) 0.2, -1 = P(left) 0.8.
B = nan(n,1);
B(abs(pLeft-.2)<1e-6)=1;
B(abs(pLeft-.5)<1e-6)=0;
B(abs(pLeft-.8)<1e-6)=-1;

% Reconstruct full-fit regressors with exactly the history conventions
% used by comparePsychometricModelsCV.
q = full.q;
k = full.block.bestK;
Bshift = nan(n,1);
for s = 1:nSessions
    ii = (starts(s):ends(s))';
    lagged = max(1,min(numel(ii),(1:numel(ii))'-k));
    Bshift(ii) = B(ii(lagged));
end
Hside = history(X,q,full.side.bestTau,starts,ends,false);
Hcontrast = history(X,q,full.contrast.bestTau,starts,ends,true);
assert(isfield(full,'drift')&&isfield(full.drift,'biasPredict'), ...
    'results.full.drift.biasPredict is required.');
Hdrift = double(full.drift.biasPredict(:));
assert(numel(Hdrift)==n,'Full drift prediction must be N-by-1.');

if strcmpi(string(opt.quantity),'bias')
    bB = full.block.beta(:);
    bS = full.side.beta(:);
    bC = full.contrast.beta(:);
    b0 = full.baseline.beta(1);
    Y = [bB(1)+bB(3)*Bshift, ...
         bS(1)+bS(3)*Hside, ...
         bC(1)+bC(3)*Hcontrast, ...
         b0+Hdrift];
    yLabel = 'Bias (log-odds)';
else
    Y = [Bshift,Hside,Hcontrast,Hdrift];
    yLabel = 'Bias regressor';
end

% Identify ONLY direct switches between the two biased block types.
% Each new block starts at t=0 and cannot cross a session boundary.
switches = cell(nSessions,2);
transitionIndices = cell(nSessions,2);
nTransitions = zeros(nSessions,2);
traces = nan(numel(rel),nSessions,4,2);
for s = 1:nSessions
    ii = (starts(s):ends(s))';
    b = B(ii);
    t02to08 = ii(find(b(1:end-1)==1 & b(2:end)==-1)+1);
    t08to02 = ii(find(b(1:end-1)==-1 & b(2:end)==1)+1);
    switches{s,1}=t02to08;
    switches{s,2}=t08to02;

    for d = 1:2
        events = switches{s,d};
        transitionIndices{s,d}=events;
        nTransitions(s,d)=numel(events);
        if isempty(events), continue; end

        for m = 1:4
            epochs = nan(numel(rel),numel(events));
            for e = 1:numel(events)
                idx = events(e)+rel;
                inside = idx>=starts(s) & idx<=ends(s);
                epochs(inside,e) = Y(idx(inside),m);
            end
            % Each session contributes one mean trace per direction/model.
            traces(:,s,m,d) = mean(epochs,2,'omitnan');
        end
    end
end

%Average traces over groups of consecutive sessions

if isempty(opt.sessionEdges)
    sessionGroups = num2cell((1:nSessions)');
    groupLabels = labels;
else
    edges = double(opt.sessionEdges(:)');

    assert(numel(edges)>=2 && all(diff(edges)>0), ...
        'sessionEdges must contain at least two increasing edges.');

    assert(edges(1)>=1 && edges(end)<=nSessions+1, ...
        'sessionEdges must lie between 1 and nSessions+1.');

    sessionGroups = cell(numel(edges)-1,1);
    groupLabels = strings(numel(edges)-1,1);

    for g = 1:numel(sessionGroups)
        sessionGroups{g} = edges(g):(edges(g+1)-1);

        groupLabels(g) = sprintf('Sessions %d-%d', ...
            edges(g),edges(g+1)-1);
    end
end

nGroups = numel(sessionGroups);

groupTraces = nan(numel(rel),nGroups,4,2);
groupNTransitions = zeros(nGroups,2);
groupCenters = nan(nGroups,1);

for g = 1:nGroups
    ss = sessionGroups{g};

    groupCenters(g) = mean(ss);
    groupNTransitions(g,:) = sum(nTransitions(ss,:),1);

    for m = 1:4
        for d = 1:2
            groupTraces(:,g,m,d) = ...
                mean(traces(:,ss,m,d),2,'omitnan');
        end
    end
end

% Chronological diverging palette: early sessions blue, late sessions red.
% A mid-tone grey (rather than white) keeps middle sessions visible on white.
% Session order follows their first occurrence in trialsT.
anchorPos = [0 0.5 1];
anchorRGB = [
    0.13 0.35 0.65;
    0.62 0.62 0.66;
    0.70 0.14 0.23];
if nSessions==1
    groupColors = anchorRGB(2,:);
else
    groupColors = interp1(anchorPos,anchorRGB, ...
        (groupCenters-1)/(nSessions-1),'linear');
end
colors = groupColors;
modelNames = {'Delayed block','Stimulus side', ...
              'Contrast weighted','Drifting bias'};
directionNames = {'P(left): 0.2 \rightarrow 0.8', ...
                  'P(left): 0.8 \rightarrow 0.2'};

fig = figure('Color','w','Position',[100 70 1150 850]);
tl = tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');
ax = gobjects(4,2);
for m = 1:4
    for d = 1:2
        ax(m,d) = nexttile(tl,(m-1)*2+d);
        hold(ax(m,d),'on');
        for g = 1:nGroups
            v = groupTraces(:,g,m,d);
            if nnz(isfinite(v)) < opt.minTrials || groupNTransitions(g,d)==0
                continue
            end
            plot(ax(m,d),rel,v,'-', ...
                'Color',colors(g,:), ...
                'LineWidth',1.15, ...
                'DisplayName',char(groupLabels(g)));
        end
        xline(ax(m,d),0,':','Color',[.35 .35 .35], ...
            'LineWidth',1.0,'HandleVisibility','off');
        xlim(ax(m,d),[rel(1) rel(end)]);
        ax(m,d).TickDir='out';
        ax(m,d).FontSize=9;
        box(ax(m,d),'off');
        if m==1
            title(ax(m,d),directionNames{d},'FontWeight','bold');
        end
        if d==1
            ylabel(ax(m,d),sprintf('%s\n%s',modelNames{m},yLabel));
        end
        if m==4
            xlabel(ax(m,d),'Trials relative to block transition');
        else
            ax(m,d).XTickLabel=[];
        end
    end
end

% Comparable y-limits across the two transition directions for each model.
for m = 1:4
    v = reshape(groupTraces(:,:,m,:),[],1);
    v = v(isfinite(v));
    if ~isempty(v)
        lo=min(v); hi=max(v);
        pad=max((hi-lo)*.08,1e-3);
        if hi==lo, pad=max(abs(hi)*.08,.1); end
        ylim(ax(m,1),[lo-pad hi+pad]);
        ylim(ax(m,2),[lo-pad hi+pad]);
    end
end

% Dummy handles ensure legend includes sessions even when a session lacks
% transitions of the particular type in the first subplot.
legendHandles = gobjects(nGroups,1);
for g = 1:nGroups
    legendHandles(g) = plot(ax(1,2),nan,nan,'-', ...
        'Color',colors(g,:), ...
        'LineWidth',1.6);
end
lgd = legend(ax(1,2),legendHandles,cellstr(groupLabels), ...
    'Location','eastoutside', ...
    'Box','off','FontSize',8);
lgd.Title.String = 'Session order (early to late)';

title(tl,sprintf('Full-fit bias around block transitions (%s)', ...
    lower(char(opt.quantity))));

out.relativeTrials = rel;
out.sessionLabels = labels;
out.colors = colors;
out.traces = traces;
out.nTransitions = nTransitions;
out.transitionIndices = transitionIndices;
out.modelNames = modelNames;
out.directionNames = directionNames;
out.biasPerTrial = Y;
out.groupTraces = groupTraces;
out.groupLabels = groupLabels;
out.sessionGroups = sessionGroups;
out.groupColors = groupColors;
out.groupNTransitions = groupNTransitions;
end

function H = history(X,q,tau,starts,ends,weighted)
% Same recursive update as historyPriors in comparePsychometricModelsCV.
H=zeros(size(X));
a=exp(-1/tau);
for s=1:numel(starts)
    num=0;den=0;
    for t=starts(s):ends(s)
        if den>0, H(t)=num/den; end
        x=X(t);
        if ~isfinite(x)||x==0,continue;end
        if weighted,w=abs(x)^q;else,w=1;end
        num=a*num+w*sign(x);
        den=a*den+w;
    end
end
end
