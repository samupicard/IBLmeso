function [pvals,Us] = getMannWhitneyU(allResps, trialCnd)

% computes Mann-Whitney U test statistic of each neuron's avg response across two sets of trials.
%
%inputs
%   allResps [ROIs x trials] avg response of each neuron to each trial
%   trialCnd [2 x trials] logical array of trials to compare w/ each other
%   NB: computes resps in LAST set in trialCnd minus resps in FIRST set in trialCnd.

%outputs
%   Us [1, ROIs] mann-whitney U test statistic
%   pvals [1, ROIs] mann-whitney U test p-value (p = U/(n1*n2))
%
% written by Samuel Picard (Nov 2023)

manualflag = false;

if isa(trialCnd,'double')
    trialCnd = logical(trialCnd);
end


if manualflag
    [~,ranks] = sort(allResps,2,'ascend');
    R1 = sum(ranks(:,trialCnd(1,:)),2);
    R2 = sum(ranks(:,trialCnd(end,:)),2);
    n1 = sum(trialCnd(1,:));
    n2 = sum(trialCnd(end,:));
    Us = min([R1 - n1*(n1+1)/2; R2 - n2*(n2+1)/2]);
    pvals = Us / n1*n2;
else
    Us = nan(1,size(allResps,1));
    pvals = nan(1,size(allResps,1));
    for i=1:size(allResps,1)
        [pvals(i),Us(i),~] = ranksum(allResps(i,trialCnd(1,:)),allResps(i,trialCnd(end,:)));
    end
end
