function R = testStimHistoryTauTrend(trialsT,taus,sessionDates,nPerm)
%TESTSTIMHISTORYTAUTREND Test changes in fitted stimulus-history timescale.
%
% R = testStimHistoryTauTrend(trialsT,taus,sessionDates,nPerm)
%
% Inputs
%   trialsT       concatenated trials table with sessionID (one ID/session)
%   taus          candidate decays, in nonzero-contrast trials (default 2:2:100)
%   sessionDates  optional datetime vector, one entry per session in stable
%                 sessionID order. If omitted uses session ordinal number.
%   nPerm         number of permutations (default 9999)
%
% Requirements: classifyEngagementTrials, fitPM_nobias,
%               computeStimulusHistory from earlier code.
%
% q is fit ONCE from pooled eligible trials with a contrast-only model,
% then fixed in each session and for both stimulus-history variants.
% q is estimated descriptively (not from held-out trials): this is an
% exploratory *across-session trend* test, not a held-out model comparison.
%
% Returns tables R.perSession and R.trend, history-model DE profiles,
% and plots. Permutations test against random ordering of whole sessions.

if nargin<2 || isempty(taus), taus = 1:1:50; end
if nargin<3, sessionDates = []; end
if nargin<4 || isempty(nPerm), nPerm = 999; end
taus = unique(double(taus(:)'));
assert(all(isfinite(taus)&taus>0) && numel(taus)>=2, ...
    'Supply at least two positive finite tau candidates.');
assert(ismember('sessionID',trialsT.Properties.VariableNames), ...
    'trialsT must contain sessionID (load as a superSession).');

trialsT = classifyEngagementTrials(trialsT);
if ismember('contrastDiff',trialsT.Properties.VariableNames)
    Xall = double(trialsT.contrastDiff(:));
else
    L = double(trialsT.contrastLeft(:));
    Rr = double(trialsT.contrastRight(:));
    L(isnan(L)) = 0; Rr(isnan(Rr)) = 0;
    Xall = Rr-L;
end
choice = -double(trialsT.choice(:));
valid = ~trialsT.noGoTrials & ~trialsT.disengagedTrials & ...
    ismember(choice,[-1 1]) & isfinite(Xall) & ...
    ismember(round(double(trialsT.probabilityLeft(:)),1),[0.2,0.8]);

% Fix q once from subject-level contrast-only model.
[coefQ,~,~] = fitPM_nobias(double(choice(valid)==1),Xall(valid),'logit');
q = coefQ.q;

[id,~,grp] = unique(trialsT.sessionID,'stable');
nS = numel(id);
if ~isempty(sessionDates)
    assert(isdatetime(sessionDates) && numel(sessionDates)==nS, ...
        'sessionDates must be a datetime vector with one date per session.');
    sessionDates = sessionDates(:);
    x = days(sessionDates-sessionDates(1));
else
    x = (1:nS)';
end
x = double(x(:));

bestTau = nan(nS,2); bestDb = nan(nS,2);
bestDE = nan(nS,2); nTrials = zeros(nS,1);
profileDE = nan(nS,numel(taus),2);
modes = {'side','contrast'};
for s=1:nS
    ids = find(grp==s);
    use = valid(ids);
    nTrials(s) = nnz(use);
    y = double(choice(ids(use))==1);
    fX = sign(Xall(ids(use))).*abs(Xall(ids(use))).^q;
    if nTrials(s)<30 || numel(unique(y))<2
        warning('Skipping session %d: insufficient valid choices.',s);
        continue
    end
    [~,devBase] = glmfit(fX,y,'binomial','link','logit');
    pNull = min(max(mean(y),eps),1-eps);
    nullDev = -2*sum(y.*log(pNull)+(1-y).*log(1-pNull));
    for m=1:2
        for j=1:numel(taus)
            % Compute history for the WHOLE session, not just included trials.
            % The history helper assumes one session when sessionID omitted.
            % Contrast model uses the same signed power transform as GLM.
            if m==1
                historyInput = Xall(ids);
                H = computeStimulusHistory(historyInput,taus(j),[], 'side');
            else
                historyInput = sign(Xall(ids)).*abs(Xall(ids)).^q;
                H = computeStimulusHistory(historyInput,taus(j),[], 'contrast');
            end
            try
                [beta,dev] = glmfit([fX,H(use)],y, ...
                    'binomial','link','logit');
                profileDE(s,j,m) = (devBase-dev)/nullDev;
                db(j) = beta(3); %#ok<AGROW>
            catch
                profileDE(s,j,m) = NaN;
                db(j) = NaN; %#ok<AGROW>
            end
        end
        [v,jBest] = max(profileDE(s,:,m));
        if isfinite(v)
            bestDE(s,m) = v;
            bestTau(s,m) = taus(jBest);
            bestDb(s,m) = db(jBest);
        end
        clear db
    end
end

R = struct();
R.q = q; R.taus = taus; R.modes = modes;
R.profileDE = profileDE;
R.perSession = table(id(:),nTrials,bestTau(:,1),bestTau(:,2), ...
    bestDE(:,1),bestDE(:,2),bestDb(:,1),bestDb(:,2), ...
    'VariableNames',{'sessionID','nTrials','tauSide','tauContrast', ...
    'deSide','deContrast','dbSide','dbContrast'});
if ~isempty(sessionDates), R.perSession.date = sessionDates; end
R.perSession.sideAtBoundary = bestTau(:,1)==taus(1) | bestTau(:,1)==taus(end);
R.perSession.contrastAtBoundary = bestTau(:,2)==taus(1) | bestTau(:,2)==taus(end);

slope = nan(2,1); p = nan(2,1); ratio = nan(2,1);
for m=1:2
    yy = log(bestTau(:,m));
    good = isfinite(yy) & isfinite(x);
    if nnz(good)<4 || numel(unique(x(good)))<3
        continue
    end
    xx = x(good); z = yy(good);
    xx0 = xx-mean(xx);
    slope(m) = sum(xx0.*(z-mean(z)))/sum(xx0.^2);
    permSlope = zeros(nPerm,1);
    for j=1:nPerm
        zp = z(randperm(numel(z)));
        permSlope(j) = sum(xx0.*(zp-mean(zp)))/sum(xx0.^2);
    end
    p(m) = (1+sum(abs(permSlope)>=abs(slope(m))))/(nPerm+1);
    ratio(m) = exp(slope(m)*(max(xx)-min(xx)));
end
R.trend = table(string(modes(:)),slope,p,ratio, ...
    'VariableNames',{'model','slopeLogTauPerX','pPermutation','tauEndOverStart'});
R.x = x;
R.xUnit = 'session number';
if ~isempty(sessionDates), R.xUnit = 'days'; end

figure('Color','w');
tiledlayout(2,1,'TileSpacing','compact');
for m=1:2
    nexttile; hold on;
    yv = bestTau(:,m);
    plot(x,yv,'ko-','LineWidth',1);
    hit = yv==taus(1) | yv==taus(end);
    plot(x(hit),yv(hit),'ro','MarkerFaceColor','r');
    good = isfinite(yv) & isfinite(x);
    if nnz(good)>=4
        fitCoef = polyfit(x(good),log(yv(good)),1);
        plot(x(good),exp(polyval(fitCoef,x(good))), 'b--','LineWidth',1.5);
    end
    set(gca,'YScale','log');
    ylabel(sprintf('%s: best \\tau',modes{m}));
    title(sprintf('%s history: permutation p = %.3g',modes{m},p(m)));
    grid on;
end
xlabel(sprintf('Recording %s',R.xUnit));

figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact');
for m=1:2
    nexttile;
    imagesc(taus,1:nS,profileDE(:,:,m));
    axis xy; colorbar;
    xlabel('Decay constant \\tau'); ylabel('Session number');
    title(sprintf('%s: incremental DE',modes{m}));
end
end
