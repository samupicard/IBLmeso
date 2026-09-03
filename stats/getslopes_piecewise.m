function [slopes_left,slopes_right] = getslopes_piecewise(allResps, trialCnd, conds, breakpoints)

% computes slopes of a linear fit of each neuron's avg response across the trial conditions specified in trialCmp
% if no breakpoint is provided, assume we only want one slope (fitted
% across all points). Output slopes_left = -slopes_right.
% if one or more breakpoints are provided, fit piecewise linear function using linear
% least squares, with continuous breakpoints. Output slopes_left from first
% breakpoint to left-most conds, output slopes_right from last
% breakpoint to right-most conds.
%
%inputs
%   allResps [nROIs, nTrials] avg response of each neuron to each trial
%   trialCnd [nCond, nTrials] logical array of condition per trial
%   conds [1, nCond] x-values of conditions (e.g. stimulus contrasts)
%   breakpoints [1, nBreakpoints] x-values of breakpoints
%
%outputs
%   slopes [nCond x neurons]
%
%
% written by Samuel Picard (Nov 2023)
% Feb 2024: deal with nans

if isa(trialCnd,'double')
    trialCnd = logical(trialCnd);
end

if nargin<3
    conds = -floor(size(trialCnd,1)/2):1:floor(size(trialCnd,1)/2);
    breakpoints = [];
elseif nargin<4
    breakpoints = [];
end

%take out nans
nantrials = isnan(allResps(1,:));
if any(nantrials)
    %warning(sprintf('Removing %d nans from trials table...',sum(nantrials)))
    allResps = allResps(:,~nantrials);
    trialCnd = trialCnd(:,~nantrials);
end

x = nan(size(allResps,2),1);
for iCond = 1:length(conds)
    x(trialCnd(iCond,:))= conds(iCond);
end
XI = [min(conds),breakpoints,max(conds)];
YI = nan(size(allResps,1),size(XI,2));

for iROI = 1:size(allResps,1)
    y = allResps(iROI,:)';
    if size(XI,2)>2
        YI(iROI,:) = lsq_lut_piecewise(x, y, XI);
    else
        YI(iROI,:) = polyfit(x,y,1);
    end
end

%figure;
%plot(x,allResps(1,:)','.',XI,YI(1,:),'+-')

slopes_left = (YI(:,1)-YI(:,2))./(XI(:,2)-XI(:,1)); %flip the sign of this slope
slopes_right = (YI(:,end)-YI(:,end-1))./(XI(:,end)-XI(:,end-1));
