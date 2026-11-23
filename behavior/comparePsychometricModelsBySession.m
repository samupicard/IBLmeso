
function out = comparePsychometricModelsBySession(trialsT,results,varargin)
%COMPAREPSYCHOMETRICMODELSBYSESSION Compare paired out-of-fold LL by session.
%
% out = comparePsychometricModelsBySession(trialsT,results, ...
%     'nBootstrap',5000,'plot',true)
%
% Accepts both older 4-model and newer 5-model CV results. Positive
% delta LL means model A predicts choices better than model B.
% Confidence intervals resample whole sessions, not individual trials.
% Blocked CV may still train on other trials from the held-out session.

ip=inputParser;
ip.addParameter('nBootstrap',5000,@(x)isscalar(x)&&x>=0&&x==round(x));
ip.addParameter('plot',true,@(x)islogical(x)&&isscalar(x));
ip.parse(varargin{:});
opt=ip.Results;

n=height(trialsT);
p=results.cvPrediction;
nModels=size(p,2);
assert(size(p,1)==n && ismember(nModels,[4 5]), ...
    'Expected 4 or 5 CV model predictions matching trialsT.');
valid=results.validTrials(:);
assert(numel(valid)==n,'validTrials has incorrect size.');

if ismember('sessionID',trialsT.Properties.VariableNames)
    ids=string(trialsT.sessionID);
elseif ismember('sessionNr',trialsT.Properties.VariableNames)
    ids=string(trialsT.sessionNr);
else
    error('No sessionID/sessionNr found: superSession must retain session identity.');
end
ids=ids(:);
[sessionNames,~,sessionCode]=unique(ids,'stable');

choice=-double(trialsT.choice(:));
y=double(choice==1);
p=max(1e-12,min(1-1e-12,p));
llTrial=y.*log(p)+(1-y).*log1p(-p);
assert(all(isfinite(llTrial(valid,:)),'all'), ...
    'Missing or nonfinite out-of-fold predictions on valid trials.');

modelNames={'baseline','block','side','contrast'};
if nModels==5
    modelNames{end+1}='drift';
end
nGroups=numel(sessionNames);
counts=accumarray(sessionCode(valid),1,[nGroups 1],@sum,0);
useGroups=find(counts>0);
N=counts(useGroups);
nSessions=numel(N);
assert(nSessions>0,'No sessions contain eligible trials.');

LL=zeros(nSessions,nModels);
for s=1:nSessions
    rows=valid & sessionCode==useGroups(s);
    LL(s,:)=sum(llTrial(rows,:),1);
end

perSession=table(sessionNames(useGroups),N, ...
    'VariableNames',{'sessionID','nTrials'});
for m=1:nModels
    perSession.(['ll_' modelNames{m}])=LL(:,m);
end

% All unique pairs, placing each history model against the baseline,
% followed by direct history-model comparisons.
pairs=[];
for m=2:nModels
    pairs=[pairs; m 1]; %#ok<AGROW>
end
for a=2:nModels
    for b=a+1:nModels
        pairs=[pairs; a b]; %#ok<AGROW>
    end
end
nPairs=size(pairs,1);
diffLL=LL(:,pairs(:,1))-LL(:,pairs(:,2));
comparison=strings(nPairs,1);
for j=1:nPairs
    a=modelNames{pairs(j,1)};
    b=modelNames{pairs(j,2)};
    comparison(j)=string(a)+" - "+string(b);
    field=[a 'Vs' upper(b(1)) b(2:end)];
    perSession.(['deltaLL_' field])=diffLL(:,j);
    perSession.(['deltaLLPerTrial_' field])=diffLL(:,j)./N;
end

totalDeltaLL=sum(diffLL,1)';
deltaLLPerTrial=totalDeltaLL/sum(N);
nSessionsPositive=sum(diffLL>0,1)';
ciLow=nan(nPairs,1);
ciHigh=nan(nPairs,1);
if opt.nBootstrap>0 && nSessions>=2
    scores=zeros(opt.nBootstrap,nPairs);
    for i=1:opt.nBootstrap
        pick=randi(nSessions,nSessions,1);
        scores(i,:)=sum(diffLL(pick,:),1)/sum(N(pick));
    end
    ci=prctile(scores,[2.5 97.5],1);
    ciLow=ci(1,:)';
    ciHigh=ci(2,:)';
end

out.perSession=perSession;
out.summary=table(comparison,totalDeltaLL,deltaLLPerTrial, ...
    nSessionsPositive,repmat(nSessions,nPairs,1),ciLow,ciHigh, ...
    'VariableNames',{'comparison','totalDeltaLL','deltaLLPerTrial', ...
    'nSessionsPositive','nSessions','ciLow','ciHigh'});
out.modelNames=string(modelNames);
if isfield(results,'driftEvaluation') && ...
        strcmpi(results.driftEvaluation,'online')
    driftNote = ['Drift predictions use earlier held-out choices ' ...
        'for sequential updating, but not the current choice. '];
else
    driftNote = ['Drift predictions update from training choices only. '];
end
out.note=['Blocked-CV session bootstrap is descriptive: a training fold ' ...
    'can contain other trials from the same recording session. ' driftNote];

if opt.plot
    out.figure=figure('Color','w','Position',[200 200 1100 600]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    x=1:nSessions;
    nexttile; hold on;
    for m=2:nModels
        plot(x,(LL(:,m)-LL(:,1))./N,'-o', ...
            'LineWidth',1.1,'DisplayName',[modelNames{m} ' - baseline']);
    end
    yline(0,'k:','HandleVisibility','off');
    ylabel('\Delta held-out LL / trial');
    title('Each history model vs contrast-only baseline');
    legend('Location','best','Box','off');
    grid on; box off;
    xticks(x); xticklabels([]);

    nexttile; hold on;
    if nModels==5
        for m=2:4
            plot(x,(LL(:,5)-LL(:,m))./N,'-o', ...
                'LineWidth',1.1,'DisplayName',['drift - ' modelNames{m}]);
        end
        title('Drifting bias versus explicit history models');
    else
        for j=1:3
            a=[2 2 3]; b=[3 4 4];
            plot(x,(LL(:,a(j))-LL(:,b(j)))./N,'-o', ...
                'LineWidth',1.1, ...
                'DisplayName',[modelNames{a(j)} ' - ' modelNames{b(j)}]);
        end
        title('Direct comparisons between history models');
    end
    yline(0,'k:','HandleVisibility','off');
    ylabel('\Delta held-out LL / trial');
    xlabel('Session');
    legend('Location','best','Box','off');
    grid on; box off;
    xticks(x);
    if nSessions<=20
        xticklabels(perSession.sessionID);
        xtickangle(45);
    end
end
end
