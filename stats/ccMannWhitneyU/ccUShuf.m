
function [ccU, ccUp, ccUSummary] = ccUShuf(resps, trialVariable, trialCondition, shufLabels)
% Returns conditions combined Mann-Whitney U statistic (a generalization of the 'ccChoiceProb'), 
% for a set of trials by combining across conditions with a decomposed mann-whitney u-stat.
%
% resps, trialVariable, trialCondition are all vectors of the same
% length. 
%
% trialVariable should have entries that are only true and false 
%
% shufLabels is a cell array with one entry per condition that appears in trialCondition. 
%   Each cell is a matrix size nx x nshuf, each column a random permutation of integers
%   from 1:nx. First column of shufLabels should be exactly 1:nx. nx is the number of "true"
%   values for that condition, and nshuf is the number of shuffle controls desired (e.g. 1000)
%
% Samuel Picard 2024 (almost literally copied from steinmetz-et-al-2019)

n = numel(resps);
nShuf = size(shufLabels{1},2)-1;

uCond = unique(trialCondition(:));
nTotal = zeros(1,1+nShuf); dTotal = 0; 
for c = 1:numel(uCond)
    inclT = trialCondition==uCond(c);
    
    chA = trialVariable & inclT;
    nA = sum(chA);
    chB = ~trialVariable & inclT;
    nB = sum(chB);
    
    n = mannWhitneyUshuf(resps(chA), ...
        resps(chB), shufLabels{c}); 
    nTotal = nTotal+n; 
    dTotal = dTotal+nA*nB;        
end

ccU = nTotal./dTotal;

t = tiedrank(ccU); 
ccUp = t(1)/(1+nShuf);  

ccUSummary = zeros(1,6);
ccUSummary(1) = ccU(1);
ccUSummary(2) = ccUp;
ccUSummary(3) = mean(ccU(2:end)); % the mean shuffle value
ccUSummary(4) = max(t)/(1+nShuf); % max significance attained by any shuffle
ccUSummary(5) = min(t)/(1+nShuf);
ccUSummary(6) = ccU(2); % ccU of just one of the shuffles
ccUSummary(7) = t(2)/(1+nShuf); % p-value of that one shuffle
