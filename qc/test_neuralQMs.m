%quick script to test the get_neuralQMs function on IBL data

datpath = 'Y:\Subjects\SP058\2024-07-30\001\alf\FOV_01';

F = readNPY(fullfile(datpath,'mpci.ROIActivityF.npy'));
Fneu = readNPY(fullfile(datpath,'mpci.ROINeuropilActivityF.npy'));
times = readNPY(fullfile(datpath,'mpci.times.npy'));

badframes = readNPY(fullfile(datpath,'mpci.badFrames.npy'));
cellClass = readNPY(fullfile(datpath,'mpciROIs.cellClassifier.npy'));
cellType = readNPY(fullfile(datpath,'mpciROIs.mpciROITypes.npy'));

[QMs, fovQMs] = get_neuralQMs(F,Fneu,times,'badframes',badframes,'iscell',logical(cellType));

%run on 1 ROI
tic
[QMs, fovQMs] = get_neuralQMs(F(:,1),Fneu(:,1),times,'badframes',badframes);
toc

%% plot outputs

%for each metric:
%histogram for all ROIs that are cells
%histogram for all ROIs that are not cells
%vertical line for average metric (across all ROIs that are cells)

figure;

subplot(2,2,1); hold on
histogram([QMs(cellType==1).noiseLevel],logspace(0,3,30),'Normalization','probability','DisplayStyle','stairs')
histogram([QMs(cellType==0).noiseLevel],logspace(0,3,30),'Normalization','probability','DisplayStyle','stairs')
xline(fovQMs.noiseLevel,'linewidth',2)
title('noise level')
set(gca,'xscale','log')

subplot(2,2,2); hold on
histogram([QMs(cellType==1).mean],linspace(-1000,2000,30),'Normalization','probability','DisplayStyle','stairs')
histogram([QMs(cellType==0).mean],linspace(-1000,2000,30),'Normalization','probability','DisplayStyle','stairs')
xline(fovQMs.mean,'linewidth',2)
title('mean')

subplot(2,2,3); hold on
histogram([QMs(cellType==1).std],linspace(0,1000,30),'Normalization','probability','DisplayStyle','stairs')
histogram([QMs(cellType==0).std],linspace(0,1000,30),'Normalization','probability','DisplayStyle','stairs')
xline(fovQMs.std,'linewidth',2)
title('std')

subplot(2,2,4); hold on
histogram([QMs(cellType==1).skew],linspace(-2,10,30),'Normalization','probability','DisplayStyle','stairs')
histogram([QMs(cellType==0).skew],linspace(-2,10,30),'Normalization','probability','DisplayStyle','stairs')
xline(fovQMs.skew,'linewidth',2)
title('skew')