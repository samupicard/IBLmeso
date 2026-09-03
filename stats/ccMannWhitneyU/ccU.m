
function [U, t] = ccU(resps, trialVariable, trialCondition)
% Returns conditions combined Mann-Whitney U statistic (a generalization of the 'ccChoiceProb'), 
% for a set of trials by combining across conditions with a decomposed mann-whitney u-stat.
%
% resps, trialVariable, trialCondition are all vectors of the same length. 
%
% trialVariable should have entries that are only true and false 
%
% this is WITHOUT shuffle controls!
%
% Samuel Picard 2024 (simplified from steinmetz-et-al-2019)

n = numel(resps);

uCond = unique(trialCondition(:));
nTotal = 0; dTotal = 0; 
for c = 1:numel(uCond)
    inclT = trialCondition==uCond(c);
    
    chA = trialVariable & inclT;
    nA = sum(chA);
    chB = ~trialVariable & inclT;
    nB = sum(chB);
    
    n = mannWhitneyU(resps(chA), ...
        resps(chB)); 
    nTotal = nTotal+n; 
    dTotal = dTotal+nA*nB;        
end

U = nTotal./dTotal;

%t = tiedrank(U); 
