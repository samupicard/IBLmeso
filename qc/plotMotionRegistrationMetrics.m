function [driftScore,excursionScore, rTime, lowFreqPowerFraction] = plotMotionRegistrationMetrics(opsFile,varargin)
%PLOTMOTIONREGISTRATIONMETRICS Plot Suite2p registration PC diagnostics.
%
% driftScore = plotMotionRegistrationMetrics(opsFile)
%
% INPUT
%   opsFile : path to Suite2p ops .mat file, or an ops struct
%
% OPTIONAL NAME-VALUE INPUTS
%
%   'nPCsPlot'        : number of registration PCs to display
%                       (default 4)
%
%   'smoothWindowMin' : width of the smoothing window used to estimate the
%                       slow PC trend, in minutes of real recording time
%                       (default 5)
%
%   'lowFreqCutoff'    : upper frequency, in Hz, used to define the
%                        low-frequency power fraction
%                        (default 0.01)
%
%   'cropCentral'     : if true, display only a central square crop of the
%                       registration PC images
%                       (default false)
%
%   'cropSize'        : width and height, in pixels, of the central crop
%                       used when cropCentral is true
%                       (default 200)
%
%   'flipImages'      : if true, animate an additional grayscale image by
%                       alternating between the low and high PC images
%                       (default true)
%
%   'flipPeriod'      : time, in seconds, between successive low/high image
%                       flips
%                       (default 0.5)
%   'saveGif'         : if true, save the animated figure as a looping GIF
%                       (default false)
%
%   'gifPath'         : full output path for the GIF. If empty, save as
%                       motionRegistrationMetrics.gif in ops.save_path
%                       (default '')
%
% OUTPUT
% driftScore : normalized measure of slow, directional drift in each PC.
%              Defined as:
%
%                  driftScore = excursionScore * abs(rTime)
%
%              where:
%
%                  excursionScore =
%                      (95th - 5th percentile of the slow trend) ...
%                      / std(raw PC timecourse)
%
%                  rTime =
%                      correlation between the slow trend and time
%
%              Thus, the score is high when the PC shows both:
%              (1) a large slow change relative to its overall variability,
%              and
%              (2) a predominantly monotonic trend over the recording.
%
%              Non-monotonic slow changes are down-weighted because their
%              correlation with time is small.
%
% Each row shows:
%   1. ops.regPC(1,iPC,:,:)              - low-PC image
%   2. ops.regPC(2,iPC,:,:)              - high-PC image
%   3. regPC(2) - regPC(1)               - difference image
%   4. ops.tPC(:,iPC) and its slow trend
%

%% Parse inputs

p = inputParser;

addParameter(p,'nPCsPlot',4,...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p,'smoothWindowMin',5,...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p,'lowFreqCutoff',0.01,...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p,'cropCentral',false,...
    @(x) islogical(x) && isscalar(x));

addParameter(p,'cropSize',200,...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p,'flipImages',true,...
    @(x) islogical(x) && isscalar(x));

addParameter(p,'flipPeriod',0.5,...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p,'saveGif',false,...
    @(x) islogical(x) && isscalar(x));

addParameter(p,'gifPath','',...
    @(x) ischar(x) || isstring(x));

parse(p,varargin{:});

opt = p.Results;


%% Load ops

if isstruct(opsFile)

    ops = opsFile;

else

    S = load(opsFile);

    if isfield(S,'ops')
        ops = S.ops;
    else
        error('File does not contain a variable called ''ops''.');
    end

end


%% Check required fields

assert(isfield(ops,'regPC'), ...
    'ops.regPC was not found.');

assert(isfield(ops,'tPC'), ...
    'ops.tPC was not found.');

assert(isfield(ops,'nrois'), ...
    'ops.nrois was not found.');

assert(isfield(ops,'badframes'), ...
    'ops.badframes was not found.');


%% Recording timing

% Actual imaging frame rate
frameRate = 32 / double(ops.nrois);

% Number of valid imaging frames
nImagingFrames = double(sum(~ops.badframes));

% Recording duration
recordingDurationSec = nImagingFrames / frameRate;

% tPC always contains a fixed number of regularly spaced samples
nTPCsamples = double(size(ops.tPC,1));

tPCsampleRate = nTPCsamples / recordingDurationSec;

% Convert desired real-time smoothing window to tPC samples
smoothWindowSec = opt.smoothWindowMin * 60;

smoothWindow = round(smoothWindowSec * tPCsampleRate);

% Sensible bounds
smoothWindow = max(3,min(smoothWindow,nTPCsamples));

% Prefer odd window length
if mod(smoothWindow,2) == 0
    smoothWindow = smoothWindow + 1;
end

% Don't exceed signal length after forcing odd
smoothWindow = min(smoothWindow,nTPCsamples);

if mod(smoothWindow,2) == 0
    smoothWindow = smoothWindow - 1;
end


%% Compute drift metrics for all 30 PCs

nPCsScore = min(30,size(ops.tPC,2));


driftScore     = nan(1,nPCsScore);
excursionScore = nan(1,nPCsScore);
rTime          = nan(1,nPCsScore);
lowFreqPowerFraction = nan(1,nPCsScore);
pcPower        = cell(1,nPCsScore);
pcFreq         = cell(1,nPCsScore);

for iPC = 1:nPCsScore

    y = double(ops.tPC(:,iPC));

    % Slow trend
    slowTrend = smoothdata(y,'movmean',smoothWindow);

    % Slow excursion relative to total PC variability
    rawSD = std(y,'omitnan');

    slowRange = prctile(slowTrend,95) - ...
        prctile(slowTrend,5);

    if rawSD > 0

        excursionScore(iPC) = slowRange / rawSD;

        % Correlation of slow trend with time:
        % high absolute values indicate a directional / monotonic trend
        t = (1:numel(slowTrend))';

        rTime(iPC) = abs(corr( ...
            t,slowTrend, ...
            'Rows','complete'));

        % Combined drift score
        driftScore(iPC) = ...
            excursionScore(iPC) * rTime(iPC);

    end

    % Power spectra
    [pcPower{iPC},pcFreq{iPC}] = getPowerSpectrum( ...
        ops.tPC(:,iPC), ...
        tPCsampleRate);

    f = pcFreq{iPC};
    p = pcPower{iPC};

    % Ignore DC
    valid = f > 0 & isfinite(p);

    f = f(valid);
    p = p(valid);

    lowIdx = f <= opt.lowFreqCutoff;

    totalPower = trapz(f,p);

    if totalPower > 0 && any(lowIdx)
        lowFreqPower = trapz(f(lowIdx),p(lowIdx));
        lowFreqPowerFraction(iPC) = lowFreqPower / totalPower;
    end

end

%% Figure layout
%
%  merge | difference | flip | ------- timecourse -------
%

fig = figure('Color','w',...
    'Name','Motion registration metrics overview');

nPCsPlot = min(opt.nPCsPlot,nPCsScore);

tl = tiledlayout(nPCsPlot,8,...
    'TileSpacing','compact',...
    'Padding','compact');

%% Figure title from ops.save_path

savePath = strrep(ops.save_path,'\','/');

% Keep path starting from subject name
tok = regexp(savePath, ...
    'Subjects/([^/]+/\d{4}-\d{2}-\d{2}/[^/]+/alf/FOV_[^/]+)', ...
    'tokens','once');

if ~isempty(tok)
    titlePath = tok{1};
else
    titlePath = savePath;
end

title(tl, ...
    sprintf('Motion registration metrics overview | %s',titlePath), ...
    'Interpreter','none');

%% Plot each PC

for iPC = 1:nPCsPlot

    %% Images

    pcLow  = squeeze(ops.regPC(1,iPC,:,:));
    pcHigh = squeeze(ops.regPC(2,iPC,:,:));

    pcDiff = pcHigh - pcLow;


    %% Optional central crop

    if opt.cropCentral

        [H,W] = size(pcLow);

        cropSize = min([opt.cropSize,H,W]);

        row0 = floor((H-cropSize)/2) + 1;
        col0 = floor((W-cropSize)/2) + 1;

        rows = row0:(row0+cropSize-1);
        cols = col0:(col0+cropSize-1);

        pcLow  = pcLow(rows,cols);
        pcHigh = pcHigh(rows,cols);
        pcDiff = pcDiff(rows,cols);

    end


    %% Robust shared scaling of low/high images

    allVals = double([pcLow(:); pcHigh(:)]);

    commonCLim = prctile(allVals,[2 98]);

    if commonCLim(1) == commonCLim(2)
        commonCLim = [min(allVals) max(allVals)];

        if commonCLim(1) == commonCLim(2)
            commonCLim = commonCLim + [-1 1];
        end
    end


    %% Normalize both images using identical limits

    lowNorm = (double(pcLow)  - commonCLim(1)) / diff(commonCLim);
    highNorm = (double(pcHigh) - commonCLim(1)) / diff(commonCLim);

    lowNorm  = max(0,min(1,lowNorm));
    highNorm = max(0,min(1,highNorm));


    %% Merge
    %
    % low  = magenta
    % high = green
    %
    % overlap -> approximately white

    mergeImg = zeros([size(pcLow) 3]);

    mergeImg(:,:,1) = lowNorm;                % red
    mergeImg(:,:,2) = highNorm;               % green
    mergeImg(:,:,3) = lowNorm;                % blue


    % --------------------------------------------------------------
    % Merge image
    % --------------------------------------------------------------

    ax1 = nexttile;

    image(ax1,mergeImg);
    axis(ax1,'image','off');

    title(ax1,sprintf('PC %d merge',iPC));


    % --------------------------------------------------------------
    % Difference image
    % --------------------------------------------------------------

    ax2 = nexttile;

    imagesc(ax2,pcDiff);
    axis(ax2,'image','off');
    colormap(ax2,gray);

    climMax = prctile(abs(double(pcDiff(:))),98);

    if climMax > 0
        clim(ax2,[-climMax climMax]);
    end

    title(ax2,'high - low');

    % --------------------------------------------------------------
    % Flipping low/high image
    % --------------------------------------------------------------

    ax3 = nexttile;

    hFlip(iPC) = imagesc(ax3,pcLow);

    axis(ax3,'image','off');
    colormap(ax3,gray);
    clim(ax3,commonCLim);

    title(ax3,'low');


    % Store both images in the image object for timer callback
    hFlip(iPC).UserData.low  = pcLow;
    hFlip(iPC).UserData.high = pcHigh;
    hFlip(iPC).UserData.state = 1;

    % --------------------------------------------------------------
    % PC timecourse
    % --------------------------------------------------------------

    ax4 = nexttile([1 4]);

    y = double(ops.tPC(:,iPC));

    slowTrend = smoothdata(y,'movmean',smoothWindow);

    timeMin = linspace( ...
        0, ...
        recordingDurationSec/60, ...
        double(nTPCsamples))';

    plot(ax4,timeMin,y, ...
        'Color',[0.75 0.75 0.75], ...
        'LineWidth',0.5);

    hold(ax4,'on');

    plot(ax4,timeMin,slowTrend, ...
        'k', ...
        'LineWidth',1.5);

    xlim(ax4,[0 recordingDurationSec/60]);

    ylabel(ax4,'PC magnitude');

    if iPC == nPCsPlot
        xlabel(ax4,'Time (min)');
    else
        ax4.XTickLabel = [];
    end

    box(ax4,'off');

    title(ax4,sprintf( ...
        'drift = %.2f   |   excursion = %.2f   |   time r = %.2f', ...
        driftScore(iPC),excursionScore(iPC),rTime(iPC)));

    % --------------------------------------------------------------
    % Power spectrum
    % --------------------------------------------------------------

    ax5 = nexttile;

    validFreq = pcFreq{iPC} > 0;

    loglog(ax5,...
        pcFreq{iPC}(validFreq),...
        pcPower{iPC}(validFreq),...
        'k',...
        'LineWidth',1);

    box(ax5,'off');

    % Set log-scale x ticks at integer exponents, label only min and max
    xL = xlim(ax5);
    expVals = ceil(log10(xL(1))) : floor(log10(xL(2)));
    xt = 10.^expVals;
    xticks(ax5,xt);
    xtlbl = repmat({''},1,numel(xt));
    if ~isempty(xtlbl)
        xtlbl{1} = sprintf('10^{%d}',expVals(1));
        xtlbl{end} = sprintf('10^{%d}',expVals(end));
    end
    xticklabels(ax5,xtlbl);
    ax5.TickLabelInterpreter = 'tex';

    xlabel(ax5,'Frequency (Hz)');
    ylabel(ax5,'Power / Hz');

    title(ax5, sprintf('PSD | low-f power = %.2f', lowFreqPowerFraction(iPC)));

end

%% Optionally save animated GIF

if opt.saveGif

    % Default save location
    if isempty(opt.gifPath)
        gifPath = fullfile(ops.save_path, ...
            'motionRegistrationMetrics.gif');
    else
        gifPath = char(opt.gifPath);
    end

    % Make sure destination folder exists
    gifDir = fileparts(gifPath);

    if ~isempty(gifDir) && ~exist(gifDir,'dir')
        mkdir(gifDir);
    end


    %% Frame 1: all low

    for i = 1:numel(hFlip)

        ud = hFlip(i).UserData;

        hFlip(i).CData = ud.low;
        title(ancestor(hFlip(i),'axes'),'low');

    end

    drawnow;

    frame = getframe(fig);
    rgb = frame2im(frame);

    [A,map] = rgb2ind(rgb,256);

    imwrite(A,map,gifPath,'gif', ...
        'LoopCount',Inf, ...
        'DelayTime',opt.flipPeriod);


    %% Frame 2: all high

    for i = 1:numel(hFlip)

        ud = hFlip(i).UserData;

        hFlip(i).CData = ud.high;
        title(ancestor(hFlip(i),'axes'),'high');

    end

    drawnow;

    frame = getframe(fig);
    rgb = frame2im(frame);

    [A,map] = rgb2ind(rgb,256);

    imwrite(A,map,gifPath,'gif', ...
        'WriteMode','append', ...
        'DelayTime',opt.flipPeriod);


    fprintf('Saved animated GIF:\n%s\n',gifPath);

end

%% Animate low/high PC images

if opt.flipImages

    fig = ancestor(tl,'figure');

    flipTimer = timer( ...
        'ExecutionMode','fixedSpacing', ...
        'Period',opt.flipPeriod, ...
        'TimerFcn',@(~,~) flipPCImages(hFlip));

    % Store timer with figure so it stays alive
    fig.UserData.flipTimer = flipTimer;

    % Stop/delete timer when figure closes
    fig.CloseRequestFcn = @(src,evt) closeMotionFigure(src);

    start(flipTimer);

end

end

function flipPCImages(hFlip)

for i = 1:numel(hFlip)

    if ~isgraphics(hFlip(i))
        continue
    end

    ud = hFlip(i).UserData;

    ax = ancestor(hFlip(i),'axes');

    if ud.state == 1

        hFlip(i).CData = ud.high;
        ud.state = 2;

        title(ax,'high');

    else

        hFlip(i).CData = ud.low;
        ud.state = 1;

        title(ax,'low');

    end

    hFlip(i).UserData = ud;

end

drawnow limitrate

end


function closeMotionFigure(fig)

% Clean up animation timer before closing figure

if isfield(fig.UserData,'flipTimer')

    t = fig.UserData.flipTimer;

    if isvalid(t)
        stop(t);
        delete(t);
    end

end

delete(fig)

end