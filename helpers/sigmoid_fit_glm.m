%fit biased sigmoid-pair to simulated psychometric data with 2 block types
%
%we try to robustly fit two sigmoids on a full set of psychometric trial-
%data with uneven numbers of trials per block type, using the function:
% y = 1/(1+exp(-v*(X+db*B-b))) 
% where db is delta-bias, b is overall bias, v is the contrast-sensitivity
% X are contrast-per-trial, B is block-per-trial, y is response-probability

%This gives us a single delta-bias term that can be used to summarize performance
%
% Samuel Picard

%% simulate some data

rng(4,"twister") %  for reproducibility

%parameters of the dataset
cvals = [-1,-.5,-.25,-.125,-.06125,0,.06125,.125,.25,.5,1]; %contrasts 
probas = [900,100]; %nr of trials per contrast-set
blocktypes = [-1,1];
sensitivity = 30; %contrast-coefficient
bias = 0.5; %bias shared between two curves
deltabias = 0.3; %difference between two bias curves

%create fake Xs (contrast sets) and Bs (block types)
x1 = [randsample(cvals(1:ceil(length(cvals)/2)),probas(1),true),randsample(cvals(ceil(length(cvals)/2):end),probas(2),true)]; 
b1= blocktypes(1)*ones(1,sum(probas)); %left blocks
x2 = [randsample(cvals(1:ceil(length(cvals)/2)),probas(2),true),randsample(cvals(ceil(length(cvals)/2):end),probas(1),true)]; 
b2= blocktypes(2)*ones(1,sum(probas)); %right blocks
X = [x1,x2]';
B = [b1,b2]';

%shuffle these
iPerm = randperm(2*sum(probas));
X = X(iPerm);
B = B(iPerm);

if false
%make 0s and 1s from a symmetric sigmoid function with noise DOESN'T WORK
y = 0.5 + 0.5*tanh(0.5 * (sensitivity*X+deltabias*sign(B)+bias)) + 0.2*randn(size(X)); %form 1: offset & scaled hyperbolic tangent function
%y = 1/(1+exp(-v*(X-deltabias*sign(B)-bias))) + 0.2*randn(size(X)); %form 2: standard logistic function
y = double(y>0.5);

else
%alternatively make a perfect sigmoid function then draw probabilistically from it
y_probas_1 = 0.5 + 0.5*tanh(0.5 * (sensitivity*cvals-deltabias+bias)); %form 1: offset & scaled hyperbolic tangent function
y_probas_2 = 0.5 + 0.5*tanh(0.5 * (sensitivity*cvals+deltabias+bias)); 
y = nan(size(X));
for i = 1:length(X)
    if B(i)<0
        y(i) = binornd(1,y_probas_1(cvals==X(i)));
    else
        y(i) = binornd(1,y_probas_2(cvals==X(i)));
    end
end
end

% figure; 
% hold on; 
% scatter(X(B>0),y(B>0),'filled'); 
% scatter(X(B<0),y(B<0),'filled'); 


%% model fitting

% we are fitting this function: y = 1/(1+exp(-(v*X^q+db*B+b))); 
% which is a logistic function with delta-bias 'db' and overall bias 'b'
F = 'logit';
[coeff,dev,stat] = fitPM_bias(y,X,B,F);

%evaluate fit
xEval = linspace(-1,1,1000)'; %get a bunch of Xs
xEval = sign(xEval).*abs(xEval).^coeff.q; %add the fitted exponent
bEval = [-ones(1,500),ones(1,500)]'; %equal sets of blocks
bEval = bEval(randperm(length(bEval))); %shuffle this
[yhat,dylo,dyhi] = glmval(stat.beta,[xEval,bEval],F,stat); %evaluate the model betas

%% plot results

%get relevant info
phats = nan(2,length(cvals));
pci_lo = nan(2,length(cvals));
pci_hi = nan(2,length(cvals));
ns = nan(2,length(cvals));
for iBlock = 1:2
    for iContrast = 1:length(cvals)
        y_sel = y(B==blocktypes(iBlock));
        x_sel = X(B==blocktypes(iBlock));
        y_sel_c = y_sel(x_sel==cvals(iContrast));
        y_sel_c_mean = mean(y_sel_c);
        r = binornd(length(y_sel_c),y_sel_c_mean);
        [phat,pci] = binofit(r,length(y_sel_c));
        ns(iBlock,iContrast) = length(y_sel_c);
        phats(iBlock,iContrast) = phat;
        pci_lo(iBlock,iContrast) = -pci(1)+phat;
        pci_hi(iBlock,iContrast) = pci(2)-phat;
        ymean(iBlock,iContrast) = y_sel_c_mean;
    end
end

figure;
hold on;
cols = colororder;
%scatter(X(B<0),y(B<0),10,cols(1,:),'.'); 
%scatter(X(B>0),y(B>0),10,cols(2,:),'.'); 

% plot means of real data
scatter(cvals,ymean(1,:),5+200*ns(1,:)/length(X),cols(1,:),'filled')
scatter(cvals,ymean(2,:),5+200*ns(2,:)/length(X),cols(2,:),'filled')
%scatter(0.99*(cvals-0.01),ymean(1,:),5+200*ns(1,:)/length(X),cols(1,:),'filled')
%scatter(0.99*(cvals+0.01),ymean(2,:),5+200*ns(2,:)/length(X),cols(2,:),'filled')
%errorbar(cvals-0.01,phats(1,:),pci_lo(1,:),pci_hi(1,:),'LineStyle','none','Color',cols(1,:));
%errorbar(cvals+0.01,phats(2,:),pci_lo(2,:),pci_hi(2,:),'LineStyle','none','Color',cols(2,:));

% plot fit lines
%plot(xEval(bEval<0),yhat(bEval<0),'Color',cols(1,:));
%plot(xEval(bEval>0),yhat(bEval>0),'Color',cols(2,:));
boundedline([xEval(bEval<0),xEval(bEval>0)],...
    [yhat(bEval<0),yhat(bEval>0)],...
    cat(3,[dylo(bEval<0),dyhi(bEval<0)],[dylo(bEval>0),dyhi(bEval>0)]),...
    'alpha');%,...
    %'Color',cols([1,2],:));

%estimated coeffs
db = coeff.db/coeff.v;
b = coeff.b/coeff.v;

xline(-b,'k')
xline(-b+db,'Color',cols(1,:));
xline(-b-db,'Color',cols(2,:));

%real coeffs from underlying distribution
xline(-bias/sensitivity,'k--')
xline(-(bias-deltabias)/sensitivity,'--','Color',cols(1,:));
xline(-(bias+deltabias)/sensitivity,'--','Color',cols(2,:));

%add delta-bias as a title
title(sprintf('v = %.1f, b = %.3f, db = %.3f',coeff.v,coeff.b,coeff.db))

xlabel('contrastDiff')
ylabel('p(R)')
