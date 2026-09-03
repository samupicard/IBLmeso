function [coeff,dev,stat] = fitPM_nobias(y,X,F,qs)

%using 'glmfit', fit a sigmoid function to psychometric data WITHOUT block
%types
% 
%The sigmoid is a standard logistic function of the following form:
% y = 1/(1+exp(-z))
% where z = v*X^q - b
%
% weights
%   b, overall bias
%   v, contrast-sensitivity (slope of the PM curve)
%   q, exponent of contrast (adjusting the contrast)
%
% predictors
%   X, 1xN vector of contrasts (typically between -1 and 1)
%
% we are trying to predict binary responses
%   y, 1xN vector of responses (0 for left, 1 for right) 
%
% Samuel Picard

% optional params
if nargin<3
    F = 'logit';
end
if nargin<5 %add exponents
    %qs = 1; 
    qs = [0.05:0.05:1];
    %qs = logspace(-2,0,15); 
end
plotQs = false; %set to true if you want to plot the q exponent

%fit the model(s)
coeffs = [];
devs = nan(1,length(qs));
stats = {};
for iMdl = 1:length(qs)
    fX = sign(X).*abs(X).^qs(iMdl);
    [coeffs(iMdl,:), devs(iMdl), stats{iMdl}] = glmfit(fX,y,'binomial','link',F);
end

%choose the best one
[~,iBestMdl] = min(devs);
q = qs(iBestMdl);
betas = coeffs(iBestMdl,:);

%plot the deviance of each model
if plotQs
    figure;
    plot(qs,devs);
    hold on;
    xline(qs(iBestMdl),'r');
    xlabel('q');
    ylabel('dev');
end

%outputs
coeff = struct('v',betas(2),'b',betas(1),'q',q);
dev = devs(iBestMdl);
stat = stats{iBestMdl};
