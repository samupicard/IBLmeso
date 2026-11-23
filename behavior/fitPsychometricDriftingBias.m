function R = fitPsychometricDriftingBias(trialsT,varargin)
%FITPSYCHOMETRICDRIFTINGBIAS Logistic PM with Gaussian random-walk bias.
%
% R = fitPsychometricDriftingBias(trialsT,'fixedQ',0.6,'plot',true)
%
% P(right_t) = logistic(v*sign(X_t)*abs(X_t)^q + b_t)
% b_1 = 0 exactly in each session
% b_t = b_(t-1) + sigma*epsilon_t, epsilon_t ~ N(0,1)
%
% v and q are fixed from a pooled contrast-only baseline fit (or q can
% be supplied via fixedQ). sigma is estimated by maximizing a Laplace
% approximation to the marginal likelihood, integrating out b_2,...,b_T.
% The MAP bias trajectory is then recovered for each session.
%
% Requires Statistics and Machine Learning Toolbox (glmfit).
% NoGo/disengagement classifier: classifyEngagementTrials.m
%
% Name-value options:
%   fixedQ       scalar or [] (default []; estimated on pooled trials)
%   sigmaBounds  [lower upper] (default [.002 .5]) in log-odds/trial
%   sigmaGrid    positive vector for optional evidence profile (default [])
%   include50    include initial 0.5-block choices in fitting (default true)
%   plot         plot example results (default true)
%
% Outputs:
%   R.sigma          fitted per-trial drift SD
%   R.q, R.v         fixed sensory exponent and coefficient
%   R.biasMAP        N x 1 smoothed MAP bias, starts at 0 per session
%   R.pMAP           N x 1 in-sample fitted probabilities (NOT CV)
%   R.validTrials    N x 1 mask of trials used as choice observations
%   R.logEvidence    Laplace log marginal likelihood (up to constants)
%   R.evidenceProfile table if sigmaGrid supplied
%
% Note: smoothing uses future choices. Do not use R.pMAP as held-out
% predictive probabilities in model comparison.

ip = inputParser;
ip.addParameter('fixedQ',[],@(x) isempty(x)||(isscalar(x)&&isfinite(x)&&x>0));
ip.addParameter('sigmaBounds',[.002 .5],@(x) isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&all(x>0)&&x(1)<x(2));
ip.addParameter('sigmaGrid',[],@(x) isempty(x)||(isnumeric(x)&&all(x(:)>0)));
ip.addParameter('plotSigmaFactors',1,@(x) isnumeric(x)&&all(isfinite(x))&&all(x>0));
ip.addParameter('initialBiasSD',0.5,@(x) isscalar(x) && isfinite(x) && x>=0);
ip.addParameter('include50',true,@(x) islogical(x)&&isscalar(x));
ip.addParameter('plot',true,@(x) islogical(x)&&isscalar(x));
ip.parse(varargin{:});
opt = ip.Results;

n = height(trialsT);
if ~ismember('noGoTrials',trialsT.Properties.VariableNames) || ...
   ~ismember('disengagedTrials',trialsT.Properties.VariableNames)
    trialsT = classifyEngagementTrials(trialsT);
end

choice = -double(trialsT.choice(:));
if ismember('contrastDiff',trialsT.Properties.VariableNames)
    X = double(trialsT.contrastDiff(:));
else
    cr = double(trialsT.contrastRight(:));
    cl = double(trialsT.contrastLeft(:));
    cr(isnan(cr))=0; cl(isnan(cl))=0;
    X = cr-cl;
end

valid = ~trialsT.noGoTrials(:) & ~trialsT.disengagedTrials(:) & ...
    ismember(choice,[-1 1]) & isfinite(X);
if ~opt.include50
    valid = valid & abs(trialsT.probabilityLeft(:)-0.5)>1e-6;
end
assert(sum(valid)>10 && numel(unique(choice(valid)))==2, ...
    'Insufficient eligible choices for fitting.');
y = double(choice==1);

% Estimate q once from the pooled contrast-only GLM, unless supplied.
if isempty(opt.fixedQ)
    qs = .05:.05:1;
    dev = nan(numel(qs),1);
    betas = cell(numel(qs),1);
    for j=1:numel(qs)
        f = sign(X(valid)).*abs(X(valid)).^qs(j);
        [betas{j},dev(j)] = glmfit(f,y(valid),'binomial','link','logit');
    end
    [~,best] = min(dev);
    q = qs(best);
    baselineBeta = betas{best};
else
    q = opt.fixedQ;
    f = sign(X(valid)).*abs(X(valid)).^q;
    baselineBeta = glmfit(f,y(valid),'binomial','link','logit');
end
v = baselineBeta(2);
fx = zeros(n,1);
fx(isfinite(X)) = sign(X(isfinite(X))).*abs(X(isfinite(X))).^q;
etaStim = v*fx;

% Handle separate sessions even if sessionID is categorical/string.
if ismember('sessionID',trialsT.Properties.VariableNames)
    [~,~,group] = unique(trialsT.sessionID,'stable');
else
    group = ones(n,1);
end
nSessions = max(group);
ids = cell(nSessions,1);
for s=1:nSessions
    ids{s} = find(group==s);
    assert(all(diff(ids{s})==1), 'Each sessionID must occupy contiguous rows.');
end

lo = log(opt.sigmaBounds(1));
hi = log(opt.sigmaBounds(2));
objective = @(ls) -totalEvidence( ...
    exp(ls),opt.initialBiasSD,ids,etaStim,y,valid);
options = optimset('TolX',1e-3,'Display','off');
[logSigma,negEvidence] = fminbnd(objective,lo,hi,options);
sigma = exp(logSigma);

biasMAP = zeros(n,1);
pMAP = nan(n,1);
for s=1:nSessions
    ii=ids{s};
    [~,b] = sessionEvidence(sigma,opt.initialBiasSD,etaStim(ii),y(ii),valid(ii));
    biasMAP(ii)=b;
    pMAP(ii)=logisticStable(etaStim(ii)+b);
end

%% Forward filtering

biasPredict = nan(n,1);
biasPredictVar = nan(n,1);

biasFilter = nan(n,1);
biasFilterVar = nan(n,1);

pPredict = nan(n,1);

for s = 1:nSessions

    ii = ids{s};

    [muP,varP,muF,varF,pP] = filterPsychometricBias( ...
        etaStim(ii),y(ii),valid(ii), ...
        sigma,opt.initialBiasSD,21);

    biasPredict(ii) = muP;
    biasPredictVar(ii) = varP;

    biasFilter(ii) = muF;
    biasFilterVar(ii) = varF;

    pPredict(ii) = pP;
end

% Sequential predictive log-likelihood
p = min(max(pPredict(valid),eps),1-eps);

logLikPredict = sum( ...
    y(valid).*log(p) + ...
    (1-y(valid)).*log1p(-p));

R = struct();
R.sigma = sigma;
R.q = q;
R.v = v;
R.baselineBias = baselineBeta(1);
R.biasMAP = biasMAP;
R.pMAP = pMAP;
R.validTrials = valid;
R.logEvidence = -negEvidence;
R.sessionGroup = group;
R.nSessions = nSessions;
R.sigmaAtBoundary = (logSigma-lo<.01)||(hi-logSigma<.01);
if ~isempty(opt.sigmaGrid)
    g = opt.sigmaGrid(:);
    ev = arrayfun(@(s) totalEvidence(s,opt.initialBiasSD,ids,etaStim,y,valid),g);
    R.evidenceProfile = table(g,ev,'VariableNames',{'sigma','logEvidence'});
else
    R.evidenceProfile = table();
end
R.initialBiasSD = opt.initialBiasSD;
R.biasPredict = biasPredict;
R.biasPredictSD = sqrt(biasPredictVar);
R.biasFilter = biasFilter;
R.biasFilterSD = sqrt(biasFilterVar);
R.pPredict = pPredict;
R.logLikPredict = logLikPredict;


if opt.plot
    % Refit bias trajectories for alternative sigma values
    factors = unique([opt.plotSigmaFactors(:); 1])';
    sigmaExamples = sigma*factors;
    trajectories = nan(n,numel(factors));

    for j = 1:numel(factors)
        if factors(j)==1
            trajectories(:,j) = biasMAP;
        else
            for s = 1:nSessions
                ii = ids{s};
                [~,bExample] = sessionEvidence( ...
                    sigmaExamples(j),opt.initialBiasSD,etaStim(ii),y(ii),valid(ii));
                trajectories(ii,j) = bExample;
            end
        end
    end

    figure('Color','w','Position',[160 180 1100 410]);
    ax = axes;
    hold(ax,'on');
    
    % Axis limits
    yValues = trajectories(:);
    yrange = max(yValues)-min(yValues);
    yrange = max(yrange,0.5);
    yLimits = [min(yValues)-0.12*yrange, ...
               max(yValues)+0.12*yrange];

    xlim(ax,[0.5 n+0.5]);
    ylim(ax,yLimits);

    % Objective block background patches
    pLeft = trialsT.probabilityLeft(:);
    edges = [1; find(diff(pLeft)~=0)+1; n+1];

    for j = 1:numel(edges)-1
        i1 = edges(j);
        i2 = edges(j+1)-1;

        if abs(pLeft(i1)-0.8)<1e-6
            col = [0.64 0.57 0.77];  % purple
        elseif abs(pLeft(i1)-0.2)<1e-6
            col = [0.94 0.76 0.49];  % orange
        else
            continue
        end

        patch(ax,[i1-.5 i2+.5 i2+.5 i1-.5], ...
            [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
            col,'FaceAlpha',0.17,'EdgeColor','none', ...
            'HandleVisibility','off');
    end

    % Bias trajectories
    h = gobjects(numel(factors),1);
    names = cell(numel(factors),1);

    for j = 1:numel(factors)
        if factors(j)==1
            col = [0 0 0];
            style = '-';
            lw = 2;
        elseif factors(j)<1
            col = [0.2 0.45 0.8];
            style = '--';
            lw = 1.2;
        else
            col = [0.8 0.3 0.2];
            style = '--';
            lw = 1.2;
        end

        h(j) = plot(ax,1:n,trajectories(:,j), ...
            'Color',col,'LineStyle',style,'LineWidth',lw);

        names{j} = sprintf('\\sigma = %.3g (%.2g x)', ...
            sigmaExamples(j),factors(j));
    end

    plot(ax,1:n,biasFilter, ...
    'Color',[0.1 0.6 0.3], ...
    'LineWidth',1.2);

    % Session boundaries
    for s = 2:nSessions
        xline(ax,ids{s}(1)-0.5,':', ...
            'Color',[0.5 0.5 0.5], ...
            'HandleVisibility','off');
    end

    yline(ax,0,'k:','HandleVisibility','off');

    xlabel(ax,'Trial');
    ylabel(ax,'Bias (log-odds)');
    title(ax,sprintf('Best \\sigma = %.4g',sigma));

    legend(ax,h,names,'Location','best','Box','off');
    box(ax,'off');
    grid(ax,'on');
    set(ax,'TickDir','out','Layer','top');

end

end


function total = totalEvidence(sigma,sigma0,ids,eta,y,valid)

total = 0;

for s = 1:numel(ids)
    ii = ids{s};
    total = total + sessionEvidence( ...
        sigma,sigma0,eta(ii),y(ii),valid(ii));
end

end


function [logEvidence,b] = sessionEvidence(sigma,sigma0,eta,y,valid)
% Gaussian random-walk bias with uncertain initial state.
%
% b(1) ~ N(0,sigma0^2), or fixed at 0 if sigma0 == 0.
% b(t)-b(t-1) ~ N(0,sigma^2).
%
% Returns Laplace log evidence and smoothed MAP trajectory.

eta = eta(:);
y = y(:);
valid = valid(:);

T = numel(eta);
b = zeros(T,1);

if T==0
    logEvidence = 0;
    return
end

% Difference operator for random-walk increments
D = spdiags(-ones(T,1),0,T,T) + ...
    spdiags(ones(T,1),1,T,T);
D = D(1:end-1,:);

Qfull = (D'*D)/sigma^2;

if sigma0>0
    % Gaussian prior on initial bias
    Qfull(1,1) = Qfull(1,1) + 1/sigma0^2;
    free = (1:T)';
    logPriorNorm = -log(sigma0)-(T-1)*log(sigma);
else
    % Original model: initial bias fixed at zero
    free = (2:T)';
    logPriorNorm = -(T-1)*log(sigma);
end

m = numel(free);

if m==0
    logEvidence = choiceLogLik(eta,y,valid);
    return
end

Q = Qfull(free,free);
obs = double(valid(free));
z = zeros(m,1);

% Newton optimization for the MAP trajectory
for iteration = 1:40

    zeta = eta(free)+z;
    p = logisticStable(zeta);

    g = obs.*(p-y(free)) + Q*z;
    w = obs.*p.*(1-p);

    H = Q + spdiags(w,0,m,m);
    delta = H\g;

    if max(abs(delta))<1e-7
        break
    end

    oldLoss = lossForZ(z,free,eta,y,valid,Q);

    scale = 1;

    while scale>1/1024
        znew = z-scale*delta;

        if lossForZ(znew,free,eta,y,valid,Q) <= oldLoss+1e-10
            break
        end

        scale = scale/2;
    end

    z = znew;
end

b(free) = z;

% Laplace approximation to marginal likelihood
p = logisticStable(eta(free)+z);
w = obs.*p.*(1-p);

H = Q + spdiags(w,0,m,m);
Rchol = chol(H);

logdetH = 2*sum(log(diag(Rchol)));

logEvidence = choiceLogLik(eta+b,y,valid) ...
    -0.5*(z'*Q*z) ...
    +logPriorNorm ...
    -0.5*logdetH;

end


function loss = lossForZ(z,free,eta,y,valid,Q)

newEta = eta;
newEta(free) = newEta(free)+z;

loss = -choiceLogLik(newEta,y,valid) + 0.5*z'*Q*z;

end


function ll=choiceLogLik(eta,y,valid)
z=eta(valid);
yv=y(valid);
% Stable Bernoulli log likelihood = y*z - log(1+exp(z))
ll=sum(yv.*z - (max(z,0)+log1p(exp(-abs(z)))));
end

function p=logisticStable(x)
p=zeros(size(x));
pos=x>=0;
p(pos)=1./(1+exp(-x(pos)));
e=exp(x(~pos));
p(~pos)=e./(1+e);
end


function [muPred,varPred,muFilt,varFilt,pPred] = ...
    filterPsychometricBias(eta,y,valid,sigma,sigma0,nQuad)
%FILTERPSYCHOMETRICBIAS Sequential Gaussian bias filtering.
%
% Predictions use only earlier choices.
% Missing/excluded observations do not update the filter.
%
% eta    : fixed sensory contribution
% y      : binary choices, 0/1
% valid  : which choices can update the filter
% sigma  : random-walk drift SD per trial
% sigma0 : initial bias SD
% nQuad  : number of Gaussian quadrature points (default 21)

if nargin<7 || isempty(nQuad)
    nQuad = 21;
end

eta = eta(:);
y = y(:);
valid = valid(:);

T = numel(eta);

muPred = zeros(T,1);
varPred = zeros(T,1);
muFilt = zeros(T,1);
varFilt = zeros(T,1);
pPred = nan(T,1);

% Standard-normal Gauss-Hermite quadrature nodes/weights
[nodes,weights] = normalQuadrature(nQuad);

mu = 0;
variance = sigma0^2;

for t = 1:T

    % Propagate state, except at the first trial
    if t>1
        variance = variance + sigma^2;
    end

    muPred(t) = mu;
    varPred(t) = variance;

    % Integration nodes for Gaussian prior
    bNodes = mu + sqrt(max(variance,0))*nodes;

    % Logistic probability at each node
    pNodes = logisticStable(eta(t)+bNodes);

    % Predictive probability before observing choice
    pPred(t) = sum(weights.*pNodes);

    if valid(t)

        % Weight prior distribution by current observation
        if y(t)==1
            likelihood = pNodes;
        else
            likelihood = 1-pNodes;
        end

        postWeights = weights.*likelihood;
        Z = sum(postWeights);

        postWeights = postWeights/Z;

        % Moment-matched Gaussian posterior
        mu = sum(postWeights.*bNodes);
        variance = sum(postWeights.*(bNodes-mu).^2);
    end

    muFilt(t) = mu;
    varFilt(t) = variance;
end

end


function [nodes,weights] = normalQuadrature(n)
% Quadrature for expectations under N(0,1).

a = sqrt(1:n-1);

J = diag(a,1)+diag(a,-1);

[V,D] = eig(J);

[nodes,order] = sort(diag(D));
weights = V(1,order)'.^2;

nodes = nodes(:);
weights = weights(:);

end
