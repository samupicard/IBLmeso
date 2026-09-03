function meandiffs = getmeandiffs(allResps, trialCnd)

% computes mean difference of each neuron's avg response across two sets of trials.
%
%inputs
%   allResps [neurons x trials] avg response of each neuron to each trial
%   trialCnd [2 x trials] logical array of trials to compare w/ each other
%   NB: computes resps in LAST set in trialCnd minus resps in FIRST set in trialCnd.

%outputs
%   meandiffs [1, neurons]
%
%
% written by Samuel Picard (Oct 2023)

if isa(trialCnd,'double')
    trialCnd = logical(trialCnd);
end
    
meandiffs = nanmean(allResps(:,trialCnd(end,:)),2) - nanmean(allResps(:,trialCnd(1,:)),2);
