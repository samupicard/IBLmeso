%this script check for motion correction traces in different FOVs from the
%same experiment

datpath = 'Y:\Subjects\SP035\2023-02-24\001';

fns = dir(fullfile(datpath,'alf','FOV*','mpci.ROIActivityF.npy')); %needs to contain mpci data
tiff_fn = dir(fullfile(datpath,'raw_imaging_*'));

nFOVs = length(fns);

xOff = [];
yOff = [];
for iFOV = 1:nFOVs
    
    reg_data_path = fns(iFOV).folder;
    
    %load ops file for accessing motion registration info
    load(fullfile(reg_data_path,'Fall.mat'),'ops');
    
    xOff(iFOV,:) = ops.xoff;
    yOff(iFOV,:) = ops.yoff;
    
end

%% plot motion reg offsets
figure;

ax = [];

ax(1) = subplot(2,6,[1:5]);
plot(xOff');

subplot(2,6,6); hold on;
for iFOV = 2:nFOVs
    scatter(xOff(1,:),xOff(iFOV,:));
end
axis square

ax(2) = subplot(2,6,[7:11]);
plot(yOff');

subplot(2,6,12); hold on;
for iFOV = 2:nFOVs
    scatter(yOff(1,:),yOff(iFOV,:));
end
axis square

linkaxes(ax,'x');
    