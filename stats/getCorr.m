function R2 = getCorr(allResps, trialVal)

% computes R-squared of each neuron's avg response against a
% trial-by-trial value (e.g. subjective bias)
%
%inputs
%   allResps [neurons x trials] avg response of each neuron to each trial
%   trialVal [1 x trials] trial values to correlate with

%outputs
%   r [1, neurons], the correlation coefficient for each neuron
%
%
% written by Samuel Picard (Dec 2024)

%TODO check some things about the input

n = size(allResps,1);

R2 = nan(1,n);
for iN = 1:size(allResps,1)
    r=corrcoef(allResps(iN,:),trialVal,'rows','complete');
    R2(iN) = (r(1,2))^2;
end