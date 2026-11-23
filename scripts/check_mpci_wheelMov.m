%% inspect relationship between mpci activity and wheel movement
%
% This script examines how population-wide mpci activity relates to
% wheel movement during an IBL session.
%
% NB: figures + metrics can be useful for QC-ing sync issues
%
% 1. Loads deconvolved ROI activity and wheel position data for one session/FOV.
% 2. Selects putative neurons using the cell classifier and performs PCA
%    on their activity to extract the dominant population activity patterns.
% 3. Computes wheel velocity and speed, and detects movement onsets
%    following periods of rest.
% 4. Quantifies the relationship between each PC and wheel movement using:
%       - R²: linear correlation between PC activity and wheel speed.
%       - |AUC|: separation of PC activity during rest vs movement.
%       - Movement-triggered PC averages (baseline-subtracted, mean ± SEM).
% 5. Produces two figures:
%       - Time courses of wheel velocity, movement onsets and PC scores.
%       - Summary plots of wheel behaviour and each PC's relationship
%         to movement (scatter, distributions, event-triggered averages).
%
% Main parameters:
%   datpath         - Session directory
%   nPCs            - Number of principal components to analyse
%   cellThreshold   - Minimum cell classifier score
%   speedThreshold  - Wheel speed threshold for defining movement
%   minRestDuration - Required rest period preceding movement onset (s)
%   eventWindow     - Window around movement onset for averaging (s)
%
% Dependencies: IBL_loadTrialsTable, readNPY, getWheelSpeed,
%               Statistics and Machine Learning Toolbox.


%datpath = 'Y:\Subjects\SP037\2023-03-09\001'; % suspected sync issues
%datpath = 'Y:\Subjects\SP075\2025-11-11\001'; % good wheel speed decoding
datpath = 'Y:\Subjects\SP054\2024-05-02\001'; % bad wheel speed decoding
%datpath = 'Y:\Subjects\SP076\2025-08-27\001';

fovNm = 'FOV_00';
blockNm = 'task_00';

nPCs = 5; %nr of PCs to plot
cellThreshold = 0.5; %threshold for suite2p cell classifier

%parameters for wheel-movement locked analysis
speedThreshold = 0.5;
minRestDuration = 0.6;  % seconds of rest preceding movement
eventWindow = [-0.6 1.5];  % seconds relative to movement onset
baselineWindow = [-0.6 -0.2];

%% Load data

trialsT = IBL_loadTrialsTable(datpath);

fovPath = fullfile(datpath,'alf',fovNm);
taskPath = fullfile(datpath,'alf',blockNm);

frameTimes = readNPY(fullfile(fovPath,'mpci.times.npy'));
wheelPos = readNPY(fullfile(taskPath,'_ibl_wheel.position.npy'));
wheelTimes = readNPY(fullfile(taskPath,'_ibl_wheel.timestamps.npy'));

tr = readNPY(fullfile(fovPath,'mpci.ROIActivityDeconvolved.npy'));
classScore = readNPY(fullfile(fovPath,'mpciROIs.cellClassifier.npy'));

% Ensure column vectors
frameTimes = double(frameTimes(:));
wheelPos = double(wheelPos(:));
wheelTimes = double(wheelTimes(:));
classScore = double(classScore(:));

% Basic dimension checks
assert(size(tr,1) == numel(frameTimes), ...
    'Number of trace samples does not match frameTimes.');

assert(size(tr,2) == numel(classScore), ...
    'Number of ROIs does not match classifier scores.');

assert(numel(wheelPos) == numel(wheelTimes), ...
    'Wheel position and timestamps have different lengths.');

assert(all(isfinite(frameTimes)) && all(diff(frameTimes) > 0), ...
    'frameTimes must be finite and strictly increasing.');

%% Select cells and compute PCA

goodCells = isfinite(classScore) & classScore > cellThreshold;

tr_good = double(tr(:,goodCells));

% Remove neurons with non-finite values or zero variance
validNeurons = all(isfinite(tr_good),1) & ...
    std(tr_good,0,1) > 0;

tr_good = tr_good(:,validNeurons);

fprintf('Using %d neurons for PCA.\n',size(tr_good,2));

assert(size(tr_good,2) >= 2, ...
    'Not enough valid neurons for PCA.');

% PCA: observations = timepoints, variables = neurons
% score is [nFrames x nComponents]
[~, score, latent, ~, explained] = pca(tr_good);

nPCs = min(nPCs,size(score,2));

%% Wheel velocity and speed

[wheelVel, wheelSpeed, wheelTimesNew] = ...
    getWheelSpeed(wheelPos,wheelTimes);

wheelVel = double(wheelVel(:));
wheelSpeed = double(wheelSpeed(:));
wheelTimesNew = double(wheelTimesNew(:));
wheelTimes = double(wheelTimes(:));

wheelSpeedFrames = nan(size(frameTimes));

avgWindow = median(diff(frameTimes));

for f = 1:numel(frameTimes)
    idx = abs(wheelTimesNew-frameTimes(f)) <= avgWindow/2;

    if any(idx)
        wheelSpeedFrames(f) = mean(wheelSpeed(idx),'omitnan');
    end
end

fprintf('Wheel support: %.1f%% of imaging frames.\n', ...
    100*mean(isfinite(wheelSpeedFrames)));

%% Feedback timestamps

feedbackTimes = [];

if ismember('feedback_times',trialsT.Properties.VariableNames)
    feedbackTimes = double(trialsT.feedback_times(:));
    feedbackTimes = feedbackTimes(isfinite(feedbackTimes));
end

%% Detect wheel movement onsets

% Identify movement periods
isMoving = wheelSpeed > speedThreshold;

% Detect threshold crossings
onsetIdx = find(diff([false; isMoving]) == 1);

% Require a preceding rest period, without crossing recording gaps
validOnsets = false(size(onsetIdx));

for k = 1:numel(onsetIdx)

    idx = onsetIdx(k);
    t = wheelTimesNew(idx);

    % Samples in preceding rest window
    preIdx = wheelTimesNew >= t-minRestDuration & wheelTimesNew < t;

    % Check that the rest window is sufficiently sampled
    nExpected = minRestDuration / median(diff(wheelTimesNew));

    validOnsets(k) = sum(preIdx) >= 0.8*nExpected && all(~isMoving(preIdx));
end

movementOnsets = wheelTimesNew(onsetIdx(validOnsets));

fprintf('Detected %d movement onsets.\n',numel(movementOnsets));

%% Extract movement-triggered PC responses

% Time axis at imaging frame rate
frameInterval = median(diff(frameTimes));
eventTime = (eventWindow(1):frameInterval:eventWindow(2))';

nEvents = numel(movementOnsets);

% [relative time x events x PCs]
eventPC = nan(numel(eventTime),nEvents,nPCs);

for k = 1:nEvents

    queryTimes = movementOnsets(k) + eventTime;

    % Interpolate PC time courses at event-relative timestamps
    eventPC(:,k,:) = reshape( ...
        interp1(frameTimes,score(:,1:nPCs), ...
        queryTimes,'linear',NaN), ...
        numel(eventTime),1,nPCs);
end

% Baseline subtraction separately for each event and PC
baselineIdx = eventTime >= baselineWindow(1) & ...
    eventTime <= baselineWindow(2);

baseline = mean(eventPC(baselineIdx,:,:),1,'omitnan');

eventPC = eventPC - baseline;

%% Compute R^2 and movement-vs-rest |AUC|

isMovingFrame = wheelSpeedFrames > speedThreshold;

R2 = nan(nPCs,1);
AUC = nan(nPCs,1);
AUCabs = nan(nPCs,1);

for i = 1:nPCs
    valid = isfinite(score(:,i)) & isfinite(wheelSpeedFrames);

    x = wheelSpeedFrames(valid);
    y = score(valid,i);

    if numel(x) >= 3 && std(x) > 0 && std(y) > 0
        r = corr(x,y);
        R2(i) = r^2;
    end

    validClass = valid & isfinite(isMovingFrame);
    pcVals = score(validClass,i);
    moving = isMovingFrame(validClass);

    if numel(pcVals) >= 3 && numel(unique(moving)) == 2
        % AUC = P(PC_movement > PC_rest), ties count as 0.5
        valsMove = pcVals(moving);
        valsRest = pcVals(~moving);

        % Pairwise-comparison definition of AUC
        C = bsxfun(@gt, valsMove(:), valsRest(:)');
        T = bsxfun(@eq, valsMove(:), valsRest(:)');
        AUC(i) = (sum(C(:)) + 0.5*sum(T(:))) / (numel(valsMove)*numel(valsRest));

        % sign-invariant effect size in [0,1]
        AUCabs(i) = 2*abs(AUC(i) - 0.5);
    end
end

%% Movement-triggered wheel velocity and speed

dtWheel = median(diff(wheelTimesNew));
eventTimeWheel = (eventWindow(1):dtWheel:eventWindow(2))';

nEvents = numel(movementOnsets);

eventVel = nan(numel(eventTimeWheel),nEvents);
eventSpeed = nan(numel(eventTimeWheel),nEvents);

for k = 1:nEvents

    queryTimes = movementOnsets(k) + eventTimeWheel;

    eventVel(:,k) = interp1( ...
        wheelTimesNew,wheelVel,queryTimes,'linear',NaN);

    eventSpeed(:,k) = interp1( ...
        wheelTimesNew,wheelSpeed,queryTimes,'linear',NaN);

    % Prevent interpolation across gaps in wheelTimesNew
    nearestDist = abs(queryTimes - interp1( ...
        wheelTimesNew,wheelTimesNew,queryTimes,'nearest',NaN));

    invalid = nearestDist > 1.5*dtWheel;

    eventVel(invalid,k) = NaN;
    eventSpeed(invalid,k) = NaN;
end

%% Figure 1: time courses only

figure('Color','w','Position',[2000 -300 1000 150+150*nPCs]);
tl1 = tiledlayout(nPCs+1,1,'TileSpacing','compact','Padding','compact');

ax = gobjects(nPCs+1,1);

% Top row: wheel velocity + movement onsets
ax(1) = nexttile;
hold on;
plot(wheelTimesNew, wheelVel, 'k');

scatter(movementOnsets, zeros(size(movementOnsets)), ...
    12, [0.85 0.2 0.2], 'filled');

ylabel('Wheel vel');
title('Wheel velocity and detected movement onsets');

% PC time courses
for i = 1:nPCs
    ax(i+1) = nexttile;
    plot(frameTimes, score(:,i), 'k');
    ylabel(sprintf('PC%d',i));

    title(sprintf('PC%d (%.1f%% var)', i, explained(i)));
end

xlabel(tl1, 'Time (s)');
title(tl1, datpath, 'Interpreter','none');

linkaxes(ax,'x');
xlim(ax(1), [frameTimes(1)+300 frameTimes(1)+400]);

%% Figure 2: wheel and PC summaries

figure('Color','w','Position',[100 100 1500 260*(nPCs+1)]);

tl2 = tiledlayout(nPCs+1,3, ...
    'TileSpacing','compact','Padding','compact');


%Top row, column 1: wheel speed distributions

nexttile(1);
hold on;

validSpeed = wheelSpeed(isfinite(wheelSpeed));

restSpeed = validSpeed(validSpeed <= speedThreshold);
moveSpeed = validSpeed(validSpeed > speedThreshold);

% Common linear bin edges
edges = linspace(0,5,100);

% Counts normalized by TOTAL number of wheel samples
nTotal = numel(validSpeed);

countsRest = histcounts(restSpeed,edges) / nTotal;
countsMove = histcounts(moveSpeed,edges) / nTotal;

colors = colororder(gca);

histogram('BinEdges',edges,'BinCounts',countsRest, ...
    'DisplayStyle','stairs', ...
    'EdgeColor',colors(1,:), ...
    'LineWidth',1.5);

histogram('BinEdges',edges,'BinCounts',countsMove, ...
    'DisplayStyle','stairs', ...
    'EdgeColor',colors(2,:), ...
    'LineWidth',1.5);

xline(speedThreshold,'k:');

xlabel('Wheel speed');
ylabel('Fraction of all samples');
title('Wheel speed distribution');

axis square;
box off;

%Top row, columns 2-3: movement-triggered wheel responses

wheelData = {eventVel,eventSpeed};
wheelLabels = {'Wheel velocity','Wheel speed'};

for j = 1:2

    nexttile(j+1);
    hold on;

    responses = wheelData{j};

    meanResp = mean(responses,2,'omitnan');

    nValid = sum(isfinite(responses),2);
    semResp = std(responses,0,2,'omitnan') ./ sqrt(nValid);

    valid = isfinite(meanResp) & isfinite(semResp);

    x = eventTimeWheel(valid);
    y = meanResp(valid);
    e = semResp(valid);

    % Shaded SEM
    fill([x;flipud(x)], ...
         [y-e;flipud(y+e)], ...
         [0.7 0.7 0.7], ...
         'EdgeColor','none', ...
         'FaceAlpha',0.4);

    % Mean response
    plot(x,y,'k','LineWidth',1.5);

    xline(0,'r--');
    yline(0,'k:');

    xlabel('Time from movement (s)');
    ylabel(wheelLabels{j});
    title(['Movement-triggered ' lower(wheelLabels{j})]);

    xlim(eventWindow);
    axis square;
    box off;
end

for i = 1:nPCs

    %Column 1: scatter + regression
    nexttile(i*3 + 1);
    hold on;

    valid = isfinite(score(:,i)) & isfinite(wheelSpeedFrames);
    x = wheelSpeedFrames(valid);
    y = score(valid,i);

    if numel(x) >= 3 && std(x) > 0 && std(y) > 0
        scatter(x,y,6,'k','filled', ...
            'MarkerFaceAlpha',0.12, ...
            'MarkerEdgeAlpha',0.12);

        p = polyfit(x,y,1);
        xFit = linspace(min(x),max(x),100);
        yFit = polyval(p,xFit);
        plot(xFit,yFit,'r','LineWidth',1.5);
    else
        text(0.5,0.5,'Insufficient data', ...
            'Units','normalized', ...
            'HorizontalAlignment','center');
    end

    xlabel('Wheel speed');
    ylabel(sprintf('PC%d',i));
    title(sprintf('R^2 = %.3f', R2(i)));
    axis square;
    box off;

    %Column 2: distributions for rest vs movement
    nexttile(i*3 + 2);
    hold on;

    valid = isfinite(score(:,i)) & isfinite(wheelSpeedFrames);
    vals = score(valid,i);
    moving = wheelSpeedFrames(valid) > speedThreshold;

    if numel(unique(moving)) == 2
        valsRest = vals(~moving);
        valsMove = vals(moving);

        edges = linspace(prctile(vals,1), prctile(vals,99), 40);

        histogram(valsRest, edges, ...
            'Normalization','pdf', ...
            'DisplayStyle','stairs', ...
            'LineWidth',1.5);

        histogram(valsMove, edges, ...
            'Normalization','pdf', ...
            'DisplayStyle','stairs', ...
            'LineWidth',1.5);

        % Labels only in the first row, top-right corner
        if i == 1
            colors = colororder(gca);

            text(0.9,0.86,'Rest', ...
                'Units','normalized', ...
                'Color',colors(1,:), ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top');

            text(0.9,0.74,'Movement', ...
                'Units','normalized', ...
                'Color',colors(2,:), ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top');
        end

    else
        text(0.5,0.5,'Insufficient data', ...
            'Units','normalized', ...
            'HorizontalAlignment','center');
    end

    xlabel(sprintf('PC%d score',i));
    ylabel('Density');
    title(sprintf('|AUC| = %.3f', AUCabs(i)));
    box off;

    %Column 3: movement-triggered average
    nexttile(i*3 + 3);
    hold on;

    responses = eventPC(:,:,i);   % [time x events]
    meanResp = mean(responses,2,'omitnan');
    nValid = sum(isfinite(responses),2);
    semResp = std(responses,0,2,'omitnan') ./ sqrt(nValid);

    validResp = isfinite(meanResp) & isfinite(semResp);

    if any(validResp)
        x = eventTime(validResp);
        y = meanResp(validResp);
        e = semResp(validResp);

        fill([x; flipud(x)], ...
            [y-e; flipud(y+e)], ...
            [0.7 0.7 0.7], ...
            'EdgeColor','none', ...
            'FaceAlpha',0.4);

        plot(eventTime, meanResp, 'k', 'LineWidth', 1.5);
        xline(0,'r--');
        yline(0,'k:');

        % robust y-limits based on full PC scale, centered on zero
        iqrVal = diff(prctile(score(:,i),[10 90]));
        if iqrVal > 0
            ylim([-1 1]*iqrVal);
        end
    else
        text(0.5,0.5,'Insufficient data', ...
            'Units','normalized', ...
            'HorizontalAlignment','center');
    end

    if i==nPCs
        xlabel('Time from movement (s)');
    end
    ylabel(sprintf('PC%d',i));
    if i==1
        title('Movement-triggered average');
    end
    xlim(eventWindow);
    axis square;
    box off;
end

title(tl2, sprintf('%s : PC summary', datpath), 'Interpreter','none');