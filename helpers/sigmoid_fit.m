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

rng(3,"twister") %  for reproducibility

%parameters of the dataset
cvals = [-1,-.5,-.25,-.125,-.06125,0,.06125,.125,.25,.5,1]; %contrasts 
probas = [500,500]; %nr of trials per contrast-set
blocktypes = [-1,1];
sensitivity = 5; %contrast-coefficient (i.e. 'threshold')
bias = 0.05; %bias shared between two curves
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
y_probas_2 = 0.5 + 0.5*tanh(0.5 * (sensitivity*cvals+deltabias+bias)); %form 1: offset & scaled hyperbolic tangent function
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

%model fitting params
mdl = '1/(1+exp(-(v*X+db*B+b)))'; %model a sigmoid function with delta-bias 'db' and overall bias 'b'
ft = fittype(mdl,'dependent',{'y'},'independent',{'X','B'},'coefficients',{'db','b','v'}); %in R2020b
%opts = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region');
%opts = fitoptions('Method','NonlinearLeastSquares');

%fit the model
%f = fit([X,B],y,ft,opts)
f = fit([X,B],y,ft)
coeffs = coeffvalues(f);
db = coeffs(1);
b = coeffs(2);

%evaluate fit
xval = linspace(-1,1,1000)';
bval = [-ones(1,500),ones(1,500)]';
bval = bval(randperm(length(bval)));
yval = feval(f,[xval,bval]);

%% plot results
figure;
hold on;
cols = colororder;
%scatter(X(B<0),y(B<0),10,cols(1,:),'.'); 
%scatter(X(B>0),y(B>0),10,cols(2,:),'.'); 
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
scatter(0.99*(cvals-0.01),ymean(1,:),5+200*ns(1,:)/length(X),cols(1,:),'filled')
scatter(0.99*(cvals+0.01),ymean(2,:),5+200*ns(2,:)/length(X),cols(2,:),'filled')
%errorbar(cvals-0.01,phats(1,:),pci_lo(1,:),pci_hi(1,:),'LineStyle','none','Color',cols(1,:));
%errorbar(cvals+0.01,phats(2,:),pci_lo(2,:),pci_hi(2,:),'LineStyle','none','Color',cols(2,:));

plot(xval(bval<0),yval(bval<0),'Color',cols(1,:));
plot(xval(bval>0),yval(bval>0),'Color',cols(2,:));

%estimated coeffs
xline(b,'k')
xline(b+db,'Color',cols(1,:));
xline(b-db,'Color',cols(2,:));

%real coeffs from underlying distribution
xline(bias,'k--')
xline(bias+deltabias,'--','Color',cols(1,:));
xline(bias-deltabias,'--','Color',cols(2,:));

%get binomial confidence intervals for each contrast


