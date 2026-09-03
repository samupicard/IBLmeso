function [slopes, R2] = getslopes(allResps, trialCmp, avg)

% computes slopes of a linear fit of each neuron's avg response across the
% trial conditions specified in trialCmp
%
%inputs
%   allResps [nROIs x nTrials] avg response of each neuron to each trial
%   trialCmp [nCond x nTrials] logical array of trials to compare w/ each other
%   avg if 'average' or 'true' it will take average values rather than

%outputs
%   slopes [nComp x nROIs]
%   R2 [nCond x nROIs]
%
%
% written by Samuel Picard (Nov 2023)
% Oct 2025: re-wrote to match usual data structuresn, added R2

if nargin<3
    avg_flag = false;
else
    if isstring(avg) || ischar(avg)
        if strcmp(avg,'average') || strcmp(avg,'mean') || strcmp(avg,'mean')
            avg_flag = true;
        else
            error('Argument (avg) not recognized')
        end
    elseif islogical(avg)
        avg_flag = avg;
    end
end

if isa(trialCmp,'double')
    trialCmp = logical(trialCmp);
end

if avg_flag
    
    vals = nan(size(allResps,1),size(trialCmp,3));
    for i = 1:size(trialCmp,1)
        vals(:,i) = nanmean(allResps(:,squeeze(trialCmp(i,:))),2);
    end
    
    x = 1:size(trialCmp,1);
    
else
    
    x_full = sum(trialCmp .* (1:size(trialCmp,1))', 1);
    iValid = x_full>0;
    x = x_full(iValid);
    vals = allResps(:,iValid);

end


%NEW: vectorized
x = x(:)'; % ensure row
X = vals;  % each row = one ROI

% Center x and data
x_mean = mean(x);
y_mean = mean(X, 2);

x_centered = x - x_mean;
y_centered = X - y_mean;

% Slope = cov(x, y) / var(x)
slopes = (y_centered * x_centered') / sum(x_centered.^2);

% Compute fitted values (optional, for R²)
yfit = slopes .* x + y_mean;

% R² = 1 - SSres / SStot
SSres = sum((X - yfit).^2, 2);
SStot = sum((X - y_mean).^2, 2);
R2 = 1 - SSres ./ SStot;


% %OLD: for-loop method
% slopes = nan(1,size(vals,1));
% R2 = nan(1,size(vals,1));
% for iROI = 1:size(vals,1)
%     
%     p = polyfit(x,vals(iROI,:),1);
%     yfit = polyval(p, x);
%     y = vals(iROI, :);
%     R2(iROI) = 1 - sum((y - yfit).^2) / sum((y - mean(y)).^2);
%     slopes(iROI) = p(1);
%     
%     %this is waaaay too slow
%     % mdl = fitlm(x, vals(iROI, :));
%     % slopes(iROI) = mdl.Coefficients.Estimate(2);
%     % R2(iROI) = mdl.Rsquared.Ordinary;
%     % RMSE(iROI) = mdl.RMSE;
% end
