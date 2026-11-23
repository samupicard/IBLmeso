function R = psychometricDriftCV(eta,y,train,starts,ends,varargin)
%PSYCHOMETRICDRIFTCV Train drift on training choices; predict sequentially.
%
% R = psychometricDriftCV(eta,y,train,starts,ends, ...
%       'sigmaBounds',[.002 .5], 'initialBiasSD',.5, ...
%       'updateMask',train)
%
% P(y_t=1) = logistic(eta_t + b_t)
% b_(first of session) ~ N(0, initialBiasSD^2)
% b_t = b_(t-1) + N(0,sigma^2)
%
% Fits sigma using a Laplace approximation to training-set marginal
% likelihood. sigma is fitted from TRAIN choices only. During forward
% filtering, choices in updateMask are assimilated AFTER their one-step
% predictions; R.pPredict(t) never uses choice y(t). For online CV, pass
% updateMask = all eligible trials, including the held-out trials.
%
% eta must contain only training-derived fixed predictor coefficients.
% It may contain a global intercept; latent bias is a zero-mean deviation.
%
% Requires no separate state-space toolbox.

ip = inputParser;
ip.addParameter('sigmaBounds',[.002 .5],@(x)isnumeric(x)&&numel(x)==2&&all(x>0)&&x(1)<x(2));
ip.addParameter('initialBiasSD',.5,@(x)isscalar(x)&&isfinite(x)&&x>0);
ip.addParameter('nQuadrature',21,@(x)isscalar(x)&&x>=5&&x==round(x));
ip.addParameter('updateMask',[],@(x)islogical(x)||isnumeric(x));
ip.parse(varargin{:});
opt=ip.Results;

eta=double(eta(:)); y=double(y(:)); train=logical(train(:));

% Missing contrasts occur only on excluded trials; use zero sensory input
% there so predictive arrays stay finite throughout the session.
assert(all(isfinite(eta(train))), ...
    'Nonfinite sensory predictor on a training trial.');
eta(~isfinite(eta))=0;

n=numel(eta);
assert(numel(y)==n && numel(train)==n,'Input sizes mismatch.');

if isempty(opt.updateMask)
    updateMask = train;  % Backward-compatible strict held-out mode
else
    updateMask = logical(opt.updateMask(:));
    assert(numel(updateMask)==n,'updateMask size mismatch.');
    assert(all(updateMask(train)), ...
        'updateMask must contain every training observation.');
end

assert(all(isfinite(y(updateMask)) & ismember(y(updateMask),[0 1])), ...
    'All updating observations must have binary choices.');
assert(numel(starts)==numel(ends) && starts(1)==1 && ends(end)==n, ...
    'Invalid session boundaries.');
assert(all(ends(1:end-1)+1==starts(2:end)), ...
    'Session boundaries must partition the rows.');

% Only sessions containing training choices contribute to the evidence.
% Trailing portions with no training observations integrate out exactly;
% truncating them substantially reduces computation.
segments = cell(numel(starts),1);

for s=1:numel(starts)
    ii=(starts(s):ends(s))';
    last=find(train(ii),1,'last');

    if ~isempty(last)
        segments{s}=ii(1:last);
    end
end

logBounds=log(opt.sigmaBounds);

objective=@(ls) -totalEvidence( ...
    exp(ls),opt.initialBiasSD,segments,eta,y,train);

[logSigma,fmin] = fminbnd( ...
    objective,logBounds(1),logBounds(2), ...
    optimset('TolX',1e-3,'Display','off'));

sigma=exp(logSigma);

%% Forward filtering

[nodes,weights]=standardNormalQuadrature(opt.nQuadrature);

pPredict=nan(n,1);
muPredict=nan(n,1);
sdPredict=nan(n,1);
muFilter=nan(n,1);
sdFilter=nan(n,1);

for s=1:numel(starts)

    ii=(starts(s):ends(s))';

    mu=0;
    variance=opt.initialBiasSD^2;

    for t=ii'

        % Predict state before observing the current choice
        if t~=ii(1)
            variance=variance+sigma^2;
        end

        muPredict(t)=mu;
        sdPredict(t)=sqrt(variance);

        % Integrate choice probability over predicted bias distribution
        bNodes=mu+sqrt(variance)*nodes;
        pNodes=logisticStable(eta(t)+bNodes);

        pPredict(t)=weights'*pNodes;

        % Assimilate current choice AFTER generating its prediction
        if updateMask(t)

            if y(t)==1
                like=pNodes;
            else
                like=1-pNodes;
            end

            posterior=weights.*like;
            posterior=posterior/sum(posterior);

            mu=posterior'*bNodes;
            variance=max(0,posterior'*((bNodes-mu).^2));
        end

        muFilter(t)=mu;
        sdFilter(t)=sqrt(variance);
    end
end

%% Outputs

R.sigma=sigma;
R.logEvidence=-fmin;

R.sigmaAtBoundary=(logSigma-logBounds(1)<.02)|| ...
    (logBounds(2)-logSigma<.02);

R.pPredict=min(1-1e-12,max(1e-12,pPredict));

R.biasPredict=muPredict;
R.biasPredictSD=sdPredict;

R.biasFilter=muFilter;
R.biasFilterSD=sdFilter;

R.initialBiasSD=opt.initialBiasSD;
R.updateMask=updateMask;

end


function total=totalEvidence(sigma,sigma0,segments,eta,y,train)

total=0;

for s=1:numel(segments)
    ii=segments{s};

    if ~isempty(ii)
        total=total+sessionEvidence( ...
            sigma,sigma0,eta(ii),y(ii),train(ii));
    end
end

end


function logEvidence=sessionEvidence(sigma,sigma0,eta,y,valid)
% Laplace approximation, integrating over ALL latent bias states.

T=numel(eta);

if T==0
    logEvidence=0;
    return
end

% Construct random-walk prior precision matrix
if T==1

    Q=sparse(1,1,1/sigma0^2,1,1);

else

    D=sparse(1:T-1,1:T-1,-1,T-1,T) + ...
      sparse(1:T-1,2:T,1,T-1,T);

    Q=(D'*D)/sigma^2;
    Q(1,1)=Q(1,1)+1/sigma0^2;

end

%% MAP optimization

b=zeros(T,1);

for iter=1:50

    p=logisticStable(eta+b);

    resid=zeros(T,1);
    resid(valid)=p(valid)-y(valid);

    w=zeros(T,1);
    w(valid)=p(valid).*(1-p(valid));

    H=Q+spdiags(w,0,T,T);
    grad=Q*b+resid;

    step=H\grad;

    if max(abs(step))<1e-7
        break
    end

    oldLoss=negLogPosterior(b,Q,eta,y,valid);
    rate=1;

    while rate>1/2048

        candidate=b-rate*step;

        if negLogPosterior(candidate,Q,eta,y,valid)<=oldLoss+1e-10
            break
        end

        rate=rate/2;
    end

    b=candidate;
end

%% Laplace approximation to marginal likelihood

p=logisticStable(eta+b);

w=zeros(T,1);
w(valid)=p(valid).*(1-p(valid));

H=Q+spdiags(w,0,T,T);
U=chol(H);

logdetH=2*sum(log(diag(U)));

% 1/2 log det(Q) = -log(sigma0)-(T-1)*log(sigma).
logEvidence=choiceLogLik(eta+b,y,valid) ...
    -0.5*(b'*Q*b) ...
    -log(sigma0) ...
    -(T-1)*log(sigma) ...
    -0.5*logdetH;

end


function loss=negLogPosterior(b,Q,eta,y,valid)

loss=-choiceLogLik(eta+b,y,valid) ...
    +0.5*(b'*Q*b);

end


function ll=choiceLogLik(z,y,valid)

z=z(valid);
y=y(valid);

ll=sum(y.*z-max(z,0)-log1p(exp(-abs(z))));

end


function p=logisticStable(z)

p=zeros(size(z));

pos=z>=0;

p(pos)=1./(1+exp(-z(pos)));

e=exp(z(~pos));
p(~pos)=e./(1+e);

end


function [nodes,weights]=standardNormalQuadrature(n)
% Quadrature for expectations under N(0,1).

a=sqrt(1:n-1);
J=diag(a,1)+diag(a,-1);

[V,D]=eig(J);

[nodes,idx]=sort(diag(D));
weights=V(1,idx)'.^2;

nodes=nodes(:);
weights=weights(:);

end
