
% subj = 'SP072'; date = '2025-08-26'; sess = '001';
%subj = 'SP058'; date = '2024-06-26'; sess = '001';
subj = 'SP058'; date = '2024-07-29'; sess = '001';

%% load the trials table

% datpath = ['Y:\Subjects' filesep subj filesep date filesep sess];
%
% fprintf('Loading trials table... ');
% trialsT = IBL_loadTrialsTable(datpath,'sync','timeline');
% if isempty(trialsT)
%     fprintf('no extracted trials found, skipping this session!\n');
%     PETH_struct = [];
%     return
% end
% fprintf('Done!\n');

subj = 'SP061';

sPaths = IBL_listSessionPaths( ...
    'root',fullfile('Y:\Subjects',subj), ...
    'protocol',{'biasedChoiceWorld'}, ...
    'trials',200);

trialsT_all = IBL_loadSuperSessionTrialsTable(sPaths);

%trialsT = trialsT_all(1:2000,:);
%trialsT = trialsT_all(1:2500,:);
%trialsT = trialsT_all(trialsT_all.sessionID==3,:);

%% fit drifting bias model (simpler version of Roy et al 2021)

%trialsT = trialsT_all(trialsT_all.sessionID<=5,:);

%R = fitPsychometricDriftingBias(trialsT,...
%    'sigmaBounds',[0.002 0.5],'sigmaGrid',logspace(-3,0,50));

% R.sigma
% R.biasMAP
% R.logEvidence

% figure;
% semilogx(R.evidenceProfile.sigma, R.evidenceProfile.logEvidence,'k-');
% xline(R.sigma,'r--');
% xlabel('\sigma_b (log-odds per trial)');
% ylabel('Log marginal likelihood');
% grid on;

%% full model comparison in one subject
% results = comparePsychometricModelsCV(trialsT, ...
%     'lags',1:20, ...
%     'taus',1:1:50, ...
%     'nFolds',5,...
%     'initialBiasSD',0.5, ...
%     'sigmaBounds',[0.002 0.5]);
trialsT = trialsT_all(ismember(trialsT_all.sessionID,1:14),:);

results = comparePsychometricModelsCV(trialsT, ...
    'lags',1:20, ...
    'taus',1:1:50, ...
    'nFolds',5, ...
    'initialBiasSD',0.5, ...
    'sigmaBounds',[0.01 0.5], ...
    'driftEvaluation','online');

disp(results.cvTable);
%disp(results.folds);

% Examine paired held-out log-likelihood differences
%comparison = comparePsychometricModelsBySession(trialsT,results);

%disp(comparison.summary);

% Visualize trials and their predictions
plotPsychometricPredictions(trialsT,results,'cv',false);

%% block-transition averages of fitted bias terms

[fig,out] = plotBlockTransitionBias(trialsT,results,'sessionEdges',1:3:15);

%% confusion matrix across all held-out trials (1 for each model)

% Actual choices (after IBL sign conversion)
choice = -double(trialsT.choice(:));
actualRight = choice == 1;

valid = results.validTrials(:);

modelNames = results.cvTable.model;

fprintf('\n%-12s %8s %8s %8s %8s %8s\n', ...
    'Model','Actual R','Pred R','Wrong L','Wrong R','Accuracy');

for m = 1:size(results.cvPrediction,2)

    p = results.cvPrediction(valid,m);

    a = actualRight(valid);
    predictedRight = p >= 0.5;

    wrongL = sum(~a & predictedRight);
    wrongR = sum(a & ~predictedRight);

    fprintf('%-12s %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n', ...
        char(modelNames(m)), ...
        100*mean(a), ...
        100*mean(predictedRight), ...
        100*wrongL/numel(a), ...
        100*wrongR/numel(a), ...
        100*mean(predictedRight==a));
end


%% Plot model performance

valid = results.validTrials(:);
actualRight = -double(trialsT.choice(:))==1;

modelNames = string(results.cvTable.model);
nModels = numel(modelNames);

% Consistent model colours (colourblind-friendly)
modelColors = [
    0.70 0.70 0.70;   % Baseline - light grey
    0.00 0.45 0.70;   % Block - blue
    0.00 0.62 0.45;   % Side - green
    0.85 0.37 0.01;   % Contrast - orange
    0.58 0.40 0.74];  % Drift - purple

edges = 0:0.1:1;
nBins = numel(edges)-1;

figure('Color','w','Position',[80 200 650 330]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

%Held-out incremental deviance explained

nexttile;

% Exclude baseline (DeltaDE = 0 by definition)
modelIdx = 2:nModels;

b = bar(1:numel(modelIdx), ...
    results.cvTable.cvDeltaDE(modelIdx), ...
    'FaceColor','flat', ...
    'EdgeColor','none');

b.CData = modelColors(modelIdx,:);

xticks(1:numel(modelIdx));
xticklabels(modelNames(modelIdx));
xtickangle(30);

ylabel('Cross-validated \DeltaDE');
title('Predictive performance');

yline(0,'k-');
box off;
set(gca,'TickDir','out','FontSize',10);

%Calibration curves

nexttile;
hold on;

for m = 1:nModels

    p = results.cvPrediction(valid,m);
    y = actualRight(valid);

    [~,~,bin] = histcounts(p,edges);

    predicted = nan(nBins,1);
    observed = nan(nBins,1);

    for i = 1:nBins
        ii = bin==i;

        if nnz(ii)>=20
            predicted(i) = mean(p(ii));
            observed(i) = mean(y(ii));
        end
    end

    plot(predicted,observed,'-o', ...
        'Color',modelColors(m,:), ...
        'MarkerFaceColor',modelColors(m,:), ...
        'MarkerSize',4, ...
        'LineWidth',1.4, ...
        'DisplayName',char(modelNames(m)));
end

% Perfect calibration reference
plot([0 1],[0 1],'k--', ...
    'LineWidth',1, ...
    'HandleVisibility','off');

axis square;
xlim([0 1]);
ylim([0 1]);

xlabel('Predicted P(right)');
ylabel('Observed P(right)');
title('Choice prob. calibration');

legend('Location','southeast','Box','off');
%grid on;
box off;
set(gca,'TickDir','out','FontSize',10);


%% plot model comparison

figure('Color','w','Position',[200 200 900 350]);
tiledlayout(1,2,'TileSpacing','compact');

% Held-out incremental deviance explained
% nexttile;
% order = unique(results.cvTable.model(2:end), 'stable'); 
% bar(categorical(results.cvTable.model(2:end), order), ...
%     results.cvTable.cvDeltaDE(2:end));
% ylabel('Cross-validated \DeltaDE');
% title('Predictive performance');
% yline(0,'k-');
% box off;

% Held-out deviance explained
nexttile;
order = unique(results.cvTable.model, 'stable'); 
bar(categorical(results.cvTable.model, order), ...
    results.cvTable.cvDevExplained);
ylabel('Cross-validated DE');
title('Predictive performance');
yline(0,'k-');
box off;

% Selected temporal parameters
nexttile;
hold on;
plot(results.folds.fold,results.folds.bestK,'ko-', ...
    'LineWidth',1.5,'DisplayName','Block lag k');
plot(results.folds.fold,results.folds.bestTauSide,'o-', ...
    'LineWidth',1.5,'DisplayName','Stimulus-side \tau');
plot(results.folds.fold,results.folds.bestTauContrast,'s-', ...
    'LineWidth',1.5,'DisplayName','Contrast-weighted \tau');

xlabel('CV fold');
ylabel('Trials');
title('Selected temporal parameters');
legend('Location','best','Box','off');
xticks(1:5);
grid on;
box off;

%% plot relationship between the bias terms of each full model

params.q = results.full.q;
params.k = results.full.block.bestK;
params.tauSide = results.full.side.bestTau;
params.tauContrast = results.full.contrast.bestTau;
driftBiasFull = results.full.drift.biasPredict;

% Compare three experimentally defined regressors
[Rcorr,R] = plotBiasRegressorCorrelations(trialsT,params,driftBiasFull);

%% do time constants change consistently across sessions?

R = testStimHistoryTauTrend(trialsT,[0.5 1:40]);

disp(R.perSession);
disp(R.trend);

%% loop over each session of each subject

f = figure;
subjects = {'SP061'}; %,'SP067'};
%subjects = {'SP054','SP058','SP060','SP061','SP067','SP072','SP075','SP076','SP081'};

for iSubj = 1:length(subjects)
    subj = subjects{iSubj};
    sPaths_all = IBL_listSessionPaths('root',['Y:\Subjects' filesep subj],'protocol',{'biasedChoiceWorld'},'trials',200);
    %sPaths_all = sPaths_all(1:min([5,length(sPaths_all)]));

    bestLags = nan(1,length(sPaths_all));
    deltabias = nan(1,length(sPaths_all));
    deBlock = nan(1,length(sPaths_all));
    pBlock = nan(1,length(sPaths_all));

    for i = 1:length(sPaths_all)
        sPaths = sPaths_all(i);

        fprintf('Loading trials tables... ');
        trialsT = IBL_loadSuperSessionTrialsTable(sPaths);
        fprintf('Done!\n');

        %% fit with a range of lags

        fprintf('Fitting PM on (pseudo-) sessions... ');
        [psychoM,RTs] = getPsycho_fitPseudosessions(trialsT,1:20,99);
        bestLags(i) = psychoM.bestLag;
        deltabias(i) = psychoM.deltabias;
        deBlock(i) = psychoM.devExplainedBlock;
        pBlock(i) = psychoM.pBlock;
        fprintf('Done!\n');

        %% fit with lag = 0
        %[psychoM_0,RTs_0] = getPsycho_fit(trialsT,true,0);

        %% plot delta-bias and incremental DE as a function of block lag k
        if false
            figure;

            yyaxis left
            plot(psychoM.lagScan.k,psychoM.lagScan.deltabias, ...
                '-o','LineWidth',1.5);
            ylabel('Delta-bias');

            yyaxis right
            plot(psychoM.lagScan.k,psychoM.lagScan.devExplainedBlock, ...
                '-s','LineWidth',1.5);
            ylabel('\DeltaDE','Interpreter','tex');

            xline(psychoM.bestLag,'k--','Best lag');
            xline(0,'k-');

            %ylim([0,0.08]);

            xlabel('subjective block lag (trials)');
            %title('Block-dependent choice bias vs. lag');
            title(sPaths{1}(13:end),'Interpreter','none');
            grid on;
        end
        %% plot best lag data and fit
        if false
            plotPMfit(psychoM,'showLegend',false);
            title(sprintf('%s, best lag = %d trials',sPaths{1}(13:end),psychoM.bestLag),'Interpreter','none');

            %plotPMfit(psychoM_0,'showLegend',false);
            %title(sprintf('%s',sPaths{1}(13:end)),'Interpreter','none');
        end
    end



    %% plot best lag and deltabias against time

    % Extract dates from session paths
    dates = datetime(regexp(sPaths_all, '\d{4}-\d{2}-\d{2}', 'match', 'once'), ...
        'InputFormat', 'yyyy-MM-dd');

    figure(f);

    sig = pBlock <= 0.05;

    % % Plot deltabias
    % ax1 = subplot(3,length(subjects),iSubj);
    % plot(dates, deltabias, 'ko-','linewidth',1);
    % ylim([-0.1,1.5])
    % %xtickangle(45);
    % %xlabel('Date');
    % if iSubj==1
    % ylabel('\delta','Interpreter','tex');
    % end
    % grid on;
    % title(subj);


    % Plot deviance explained by subjective block
    ax1 = subplot(2,length(subjects),iSubj);
    %plot(dates, deBlock, 'ko-','linewidth',1);
    plot(deBlock, 'ko-','linewidth',1);
    hold on;
    plot(find(sig),deBlock(sig),'ko','MarkerFaceColor','k');
    ylim([0,.21])
    %xtickangle(45);
    %xlabel('Date');
    if iSubj==1
    ylabel('\DeltaDE','Interpreter','tex');
    end
    grid on;
    title(subj);

    % Plot best lag
    ax2 = subplot(2,length(subjects),length(subjects)+iSubj);
    %plot(dates, bestLags, 'ks-','linewidth',1);
    plot(bestLags, 'ks-','linewidth',1);
    hold on;
    plot(find(sig),bestLags(sig),'ks','MarkerFaceColor','k');
    %xtickangle(45);
    ylim([0,20])
    xlabel('Session');
    if iSubj==1
    ylabel('best lag (trials)');
    end
    grid on;

    linkaxes([ax1 ax2],'x');
    ax1.XTick = ax2.XTick;
    ax1.XTickLabel = [];


end

if false

%% Loop over each subject's supersession and take sliding windows of consecutive trials

%subjects = {'SP061'};
subjects = {'SP054','SP058','SP060','SP061','SP067','SP072','SP075','SP076','SP081'};

N = 1000;           % Valid trials per window
step = 500;        % Sliding window step
ks = 1:20;
nPseudo = 99;

plotLagScans = false;
plotPMfits = false;

f = figure;

for iSubj = 1:length(subjects)

    subj = subjects{iSubj};

    %% Load superSession

    sPaths_all = IBL_listSessionPaths( ...
        'root',['Y:\Subjects' filesep subj], ...
        'protocol',{'biasedChoiceWorld'},'trials',200);

    fprintf('%s: Loading %d sessions... ',subj,numel(sPaths_all));
    trialsT = IBL_loadSuperSessionTrialsTable(sPaths_all);
    fprintf('Done!\n');

    %% Identify valid trials for sliding windows

    trialsT = classifyEngagementTrials(trialsT);

    valid = ~trialsT.noGoTrials & ...
            ~trialsT.disengagedTrials & ...
            ismember(round(trialsT.probabilityLeft,1),[0.2 0.8]) & ...
            ismember(trialsT.choice,[-1 1]);

    if ismember('contrastDiff',trialsT.Properties.VariableNames)
        valid = valid & isfinite(trialsT.contrastDiff);
    else
        valid = valid & ...
            (isfinite(trialsT.contrastLeft) | ...
             isfinite(trialsT.contrastRight));
    end

    validIdx = find(valid);
    nValid = numel(validIdx);

    % Sliding windows in valid-trial coordinates
    starts = 1:step:(nValid-N+1);
    nWindows = numel(starts);

    fprintf('%s: %d valid trials, %d windows\n', ...
        subj,nValid,nWindows);

    %% Initialize outputs

    bestLags = nan(1,nWindows);
    deltabias = nan(1,nWindows);
    deBlock = nan(1,nWindows);
    pBlock = nan(1,nWindows);

    windowCenters = nan(1,nWindows);

    %% Loop over sliding windows

    for i = 1:nWindows

        % Indices into the original full trials table
        winIdx = validIdx(starts(i):starts(i)+N-1);

        % Window centre in valid-trial coordinates
        windowCenters(i) = starts(i)+(N-1)/2;

        % Preserve ALL trials and block structure.
        % Exclude choices outside the current window from fitting.
        winT = trialsT;
        keep = false(height(trialsT),1);
        keep(winIdx) = true;
        winT.choice(~keep) = NaN;

        fprintf('%s: Window %d/%d (valid trials %d-%d)\n', ...
            subj,i,nWindows,starts(i),starts(i)+N-1);

        %% Fit empirical and pseudo-session models

        [psychoM,RTs] = getPsycho_fitPseudosessions( ...
            winT,ks,nPseudo);

        bestLags(i) = psychoM.bestLag;
        deltabias(i) = psychoM.deltabias;
        deBlock(i) = psychoM.devExplainedBlock;
        pBlock(i) = psychoM.pBlock;

        %% Optional: fit with lag = 0

        if plotPMfits
            [psychoM_0,RTs_0] = getPsycho_fit(winT,true,0);
        end

        %% Optional: plot lag scan

        if plotLagScans
            figure;

            yyaxis left
            plot(psychoM.lagScan.k,psychoM.lagScan.deltabias, ...
                '-o','LineWidth',1.5);
            ylabel('Delta-bias');

            yyaxis right
            plot(psychoM.lagScan.k,psychoM.lagScan.devExplainedBlock, ...
                '-s','LineWidth',1.5);
            ylabel('\DeltaDE','Interpreter','tex');

            xline(psychoM.bestLag,'k--','Best lag');
            xline(0,'k-');

            %ylim([0,0.08]);

            xlabel('subjective block lag (trials)');
            title(sprintf('%s, window %d',subj,i));
            grid on;
        end

        %% Optional: plot psychometric fit

        if plotPMfits
            plotPMfit(psychoM);
            title(sprintf('%s, window %d, best lag = %d', ...
                subj,i,psychoM.bestLag));
        end

    end

    %% Plot sliding-window results

    figure(f);

    sig = pBlock <= 0.05;

    % Plot deviance explained by subjective block
    ax1 = subplot(2,length(subjects),iSubj);

    plot(windowCenters,deBlock,'ko-','linewidth',1);
    hold on;
    plot(windowCenters(sig),deBlock(sig),'ko','MarkerFaceColor','k');

    ylim([0,.21])
    if iSubj==1
        ylabel('\DeltaDE','Interpreter','tex');
    end
    grid on;
    title(subj);

    % Plot best lag
    ax2 = subplot(2,length(subjects),length(subjects)+iSubj);

    plot(windowCenters,bestLags,'ks-','linewidth',1);
    hold on;
    plot(windowCenters(sig),bestLags(sig),'ks','MarkerFaceColor','k');

    ylim([0,20])
    xlabel('Valid trial number');
    if iSubj==1
        ylabel('best lag (trials)');
    end
    grid on;

    linkaxes([ax1 ax2],'x');
    ax1.XTick = ax2.XTick;
    ax1.XTickLabel = [];

end
end