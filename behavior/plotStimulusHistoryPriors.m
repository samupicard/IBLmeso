function plotStimulusHistoryPriors(trialsT,results,trialRange)
%PLOTSTIMULUSHISTORYPRIORS Compare objective and fitted stimulus priors.
%
% plotStimulusHistoryPriors(trialsT,results)
% plotStimulusHistoryPriors(trialsT,results,trialRange)
%
% Inputs:
%   trialsT    : IBL trials table
%   results    : output of getPsycho_fitStimHistory
%   trialRange : optional trial indices to display (default: all)
%
% Shows:
%   Objective prior: 1 - 2*probabilityLeft
%   Side-only history: exponentially weighted stimulus side
%   Contrast-weighted history: exponentially weighted contrast
%
% Both history regressors are computed using their fitted bestTau.

if nargin<3 || isempty(trialRange)
    trialRange = 1:height(trialsT);
end

trialRange = trialRange(:)';

% Objective prior, positive = rightward
objectivePrior = 1 - 2*trialsT.probabilityLeft(:);

% Fitted history regressors
sidePrior = results.side.history(:);
contrastPrior = results.contrast.history(:);

% Extract requested range
x = trialRange;
objectivePrior = objectivePrior(x);
sidePrior = sidePrior(x);
contrastPrior = contrastPrior(x);

%% Plot

figure('Color','w','Position',[200 200 1000 350]);
hold on;

% Colors
colObjective = [0.25 0.25 0.25];
colSide = [0.15 0.45 0.75];
colContrast = [0.85 0.4 0.2];

% Objective block prior
stairs(x,objectivePrior, ...
    'Color',colObjective, ...
    'LineWidth',2);

% Stimulus-history priors
plot(x,sidePrior, ...
    'Color',colSide, ...
    'LineWidth',1.5);

plot(x,contrastPrior, ...
    'Color',colContrast, ...
    'LineWidth',1.5);

xlabel('Trial');
ylabel('Prior (rightward positive)');
ylim([-1 1]);

legend({ ...
    'Objective block prior', ...
    sprintf('Stimulus side (\\tau = %.1f)',results.side.bestTau), ...
    sprintf('Contrast weighted (\\tau = %.1f)',results.contrast.bestTau)}, ...
    'Location','best','Box','off');

grid on;
box off;

end
