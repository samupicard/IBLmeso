
dDir = 'Y:\Subjects\results\taskVariableTuning_SingleCells\chronic';

%% get the taskTuned statistics


%% run this script to produce the aggregate tables

IBL_densityMaps_batch;


%% 
out_st = roi_stat_cosine_maps(fullfile(dDir,'allROIs_ccu_stimOn_0to400_stimSide100.parquet'),150, 20, 10);
out_ch = roi_stat_cosine_maps(fullfile(dDir,'allROIs_ccu_choiceMovement_-400to0_choice.parquet'),150, 20, 10);
out_fb = roi_stat_cosine_maps(fullfile(dDir,'allROIs_ccu_feedback_0to400_feedbackType.parquet'),150, 20, 10);

%% plot cosine map gui

