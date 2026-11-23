function R = testStimHistoryTauChange(trialsT,taus,nBoot,minSessionsPerSide)
%TESTSTIMHISTORYTAUCHANGE Test whether stimulus-history tau changes across sessions.
%
% R = testStimHistoryTauChange(trialsT,taus,nBoot,minSessionsPerSide)
%
% Models (separately for side and contrast-weighted history):
%   H0: one common tau across sessions; each session has its own b,v,db.
%   H1: one tau before and another after an unknown change point;
%       each session still has its own b,v,db.
%
% Fits a single pooled contrast-only q, then holds q fixed. Tests the
% improvement in summed binomial deviance by parametric bootstrap under H0.
% Bootstrap simulates choices at each session's actual contrast/history
% sequence using fitted H0 probabilities, then repeats the entire scan.
%
% Requires: classifyEngagementTrials, fitPM_nobias, glmfit.
% tau is in nonzero-contrast presentations. The current stimulus is NOT
% included in its history, and histories restart at each session boundary.
%
% Defaults: taus=2:2:100, nBoot=199, minSessionsPerSide=2.
%
% R.(mode).changeP     one-sided bootstrap p-value for an early/late change
% R.(mode).commonTau   best constant tau
% R.(mode).earlyTau, lateTau, splitAfterSession
% R.(mode).profileDev  nSessions x nTaus profile deviances
% R.(mode).bestTau, ciLow, ciHigh (approx. 95% profile intervals)
%
% The test is conditional on q, chosen tau grid and model assumptions.

if nargin<2 || isempty(taus), taus = 2:2:100; end
if nargin<3 || isempty(nBoot), nBoot = 199; end
if nargin<4 || isempty(minSessionsPerSide), minSessionsPerSide = 2; end
taus = unique(double(taus(:)'));
assert(all(isfinite(taus)&taus>0) && numel(taus)>1);
assert(ismember('sessionID',trialsT.Properties.VariableNames), ...
    'Input must be a superSession trials table with sessionID.');
T = classifyEngagementTrials(trialsT);
if ismember('contrastDiff',T.Properties.VariableNames)
    X = double(T.contrastDiff(:));
else
    L = double(T.contrastLeft(:)); RR = double(T.contrastRight(:));
    L(isnan(L))=0; RR(isnan(RR))=0; X=RR-L;
end
choice = -double(T.choice(:));
valid = ~T.noGoTrials & ~T.disengagedTrials & ...
    ismember(choice,[-1 1]) & isfinite(X) & ...
    ismember(round(double(T.probabilityLeft(:)),1),[.2 .8]);
[~,~,group] = unique(T.sessionID,'stable');
nS = max(group);
assert(nS >= 2*minSessionsPerSide, 'Not enough recording sessions.');
[qc,~,~] = fitPM_nobias(double(choice(valid)==1),X(valid),'logit');
q = qc.q;
% A fixed, common q determines both contrast predictor and history weights.
fXall = sign(X).*abs(X).^q;
modeNames = {'side','contrast'};
R = struct(); R.q = q; R.taus = taus; R.nBoot = nBoot;
R.nSessions = nS;
for m=1:2
    Z = cell(nS,1); ys = cell(nS,1);
    nUsed = zeros(nS,1);
    for s=1:nS
        ids = find(group==s);
        use = valid(ids);
        ys{s} = double(choice(ids(use))==1);
        nUsed(s) = nnz(use);
        assert(nUsed(s)>=30 && numel(unique(ys{s}))==2, ...
            'Session %d has insufficient valid choice data.',s);
        fX = fXall(ids(use));
        c = X(ids);
        if m==2, c = sign(c).*abs(c).^q; end
        % Compute every candidate history over the complete session.
        % r and w are accumulated BEFORE the current stimulus is added.
        H = zeros(numel(ids),numel(taus));
        for j=1:numel(taus)
            a=exp(-1/taus(j)); num=0; den=0;
            for t=1:numel(ids)
                if den>0, H(t,j)=num/den; end
                if isfinite(c(t)) && c(t)~=0
                    if m==1, w=1; else, w=abs(c(t)); end
                    num=a*num+w*sign(c(t));
                    den=a*den+w;
                end
            end
        end
        Z{s}=cell(numel(taus),1);
        for j=1:numel(taus)
            Z{s}{j} = [fX,H(use,j)];
        end
    end

    % Observed session-wise likelihood profiles.
    [D,~] = fitProfiles(ys,Z);
    [bestTau,lo,hi] = profileIntervals(D,taus);
    obs = changeStatistic(D,minSessionsPerSide);

    % Fit generating probabilities under H0 at the common tau.
    p0 = cell(nS,1);
    for s=1:nS
        [b,~] = glmfit(Z{s}{obs.commonIndex},ys{s}, ...
            'binomial','link','logit');
        eta = [ones(nUsed(s),1),Z{s}{obs.commonIndex}]*b;
        p0{s}=1./(1+exp(-eta));
        p0{s}=min(max(p0{s},1e-9),1-1e-9);
    end

    statBoot=nan(nBoot,1);
    for b=1:nBoot
        yBoot = cell(nS,1);
        for s=1:nS
            yBoot{s}=double(rand(nUsed(s),1)<p0{s});
        end
        Db=fitProfiles(yBoot,Z);
        if all(isfinite(Db(:)))
            tmp=changeStatistic(Db,minSessionsPerSide);
            statBoot(b)=tmp.improvement;
        end
        if mod(b,25)==0 || b==nBoot
            fprintf('%s: bootstrap %d/%d\n',modeNames{m},b,nBoot);
        end
    end
    nValid=sum(isfinite(statBoot));
    assert(nValid==nBoot, ...
        'Nonfinite bootstrap fits; reduce tau range or inspect choices.');
    p=(1+sum(statBoot>=obs.improvement))/(nBoot+1);

    S = struct();
    S.profileDev=D; S.bestTau=bestTau; S.ciLow=lo; S.ciHigh=hi;
    S.nTrials=nUsed; S.commonTau=taus(obs.commonIndex);
    S.earlyTau=taus(obs.earlyIndex);
    S.lateTau=taus(obs.lateIndex);
    S.splitAfterSession=obs.split;
    S.devianceImprovement=obs.improvement;
    S.bootstrapImprovement=statBoot; S.changeP=p;
    R.(modeNames{m})=S;
end

% Show session-wise tau estimates and their profile uncertainty.
figure('Color','w');
tiledlayout(2,1,'TileSpacing','compact');
for m=1:2
    S=R.(modeNames{m});
    nexttile; hold on;
    s=(1:nS)';
    errorbar(s,S.bestTau,S.bestTau-S.ciLow,S.ciHigh-S.bestTau, ...
        'ko','MarkerFaceColor','w','CapSize',3);
    xline(S.splitAfterSession+0.5,'k:');
    yline(S.commonTau,'Color',[.5 .5 .5],'LineStyle','--');
    plot([1 S.splitAfterSession],[S.earlyTau S.earlyTau], ...
        '-','LineWidth',2);
    plot([S.splitAfterSession+1 nS],[S.lateTau S.lateTau], ...
        '-','LineWidth',2);
    ylabel(sprintf('%s: \tau',modeNames{m}));
    title(sprintf('Change-point bootstrap p = %.3g',S.changeP));
    ylim([min(taus)*.8 max(taus)*1.1]); xlim([.5 nS+.5]);
    grid on;
end
xlabel('Recording session');
end

function [D,B] = fitProfiles(ys,Z)
nS=numel(ys); nT=numel(Z{1});
D=nan(nS,nT); B=cell(nS,nT);
for s=1:nS
    for j=1:nT
        try
            [B{s,j},D(s,j)] = glmfit(Z{s}{j},ys{s}, ...
                'binomial','link','logit');
        catch
            D(s,j)=NaN;
        end
    end
end
end

function O=changeStatistic(D,minS)
% Sum session-wise deviance; shared tau under null, two taus under alt.
[nS,~]=size(D);
[commonD,cIdx]=min(sum(D,1));
changeD=Inf; split=NaN; eIdx=NaN; lIdx=NaN;
for c=minS:nS-minS
    [eD,e]=min(sum(D(1:c,:),1));
    [lD,l]=min(sum(D(c+1:end,:),1));
    if eD+lD<changeD
        changeD=eD+lD; split=c; eIdx=e; lIdx=l;
    end
end
O.commonIndex=cIdx; O.earlyIndex=eIdx; O.lateIndex=lIdx;
O.split=split; O.improvement=max(0,commonD-changeD);
end

function [best,lo,hi]=profileIntervals(D,taus)
nS=size(D,1);
best=nan(nS,1); lo=best; hi=best;
for s=1:nS
    ds=D(s,:);
    [v,i]=min(ds);
    if ~isfinite(v), continue; end
    best(s)=taus(i);
    accepted=find(ds-v<=3.841459);
    if ~isempty(accepted)
        lo(s)=taus(accepted(1)); hi(s)=taus(accepted(end));
    end
end
end
