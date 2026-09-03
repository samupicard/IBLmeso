# IBL mesoscope analysis

WIP analyses of IBL 2p-RAM data. All Matlab based.

## Requirements

- MATLAB R2022b or later (tested on R2025b)
- Required toolboxes:
  - Statistics and Machine Learning Toolbox
  - Parallel Computing Toolbox
- npy-matlab
- glmnet_matlab
- allenCCF

## Installation

Clone the repository:

```bash
git clone https://github.com/samupicard/IBLmeso.git
```

Then add the repository to your MATLAB path:

```matlab
addpath(genpath('path/to/IBLmeso'))
```

## Usage

Example 1: minimal loading

```matlab
% load paths of canonical sessions
load('canonicalSessions.mat');

%take one example session
datpath = sPaths{109}; %'Y:\Subjects\SP072\2025-08-26\001'
splitPath = split(datpath,'\');
subj = splitPath{end-2};
date = splitPath{end-1};
sess = splitPath{end};

% load neural data (deconvolved spikes) from all FOVs
Fall = IBL_loadMesoData(subj, date, sess, 'trace', 'spks', 'fast', true);

% consider all non-zero frameQC as bad
badframes = find(Fall.frameQC~=0);

% load trials table
trialsT = IBL_loadTrialsTable(datpath,'sync','timeline');

```


Example 2: compute stimulus-side PETHs in each ROI of one session

```matlab

params = struct( ...
    'activity_type',   'deconv', ...
    'trialTypeField',  'contrastDiff', ...
    'trialTypeVals',   [-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1], ...
    'trialTypeFilter', '', ...
    'condFields',      {''}, ...
    'evnt',            {'stimOn'}, ...
    'twin_all',        [-1,3], ...
    'twin_bl',         'none', ...
    'twin_ev',         {[0,0.4]}, ...
    'stat_to_use',     {'ccMeanDiff'}, ...
    'nComp',           5, ...
    'nTrialsMin',      20, ...
    'nTrialsToKeep',   false, ...
    'minTrialsPerCond',10, ...
    'cv',              '', ...
    'pthresh',         0.05, ...
    'nPseudoSessions', 199 ...
    );

PETH_struct = get_taskTunedROIs(datpath,params,'overwriteExisting', false, 'saveToFOV',false); 

```

## Data

This repository does not include raw data. 
Analysis scripts expect data to be locally accessible under 'Y:\Subjects'


## Authors

Samuel Picard