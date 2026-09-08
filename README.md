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


Example 2: For each ROI of one session, compute mean PETHs for each signed contrast, locked to stimOn, controlling for choice. Get tuning stat in time-period 0 to 400ms, for -100% v. +100%.

```matlab

params = struct( ...
    'activity_type',   'deconv', ...
    'trialTypeField',  'contrastDiff', ...
    'trialTypeVals',   [-1,-.5,-.25,-.125,-.0625,0,.0625,.125,.25,.5,1], ...
    'compareVals',     [-1, 1], ...
    'trialTypeFilter', '', ...
    'condFields',      {'choice'}, ...
    'evnt',            {'stimOn'}, ...
    'twin_all',        [-1,3], ...
    'twin_bl',         'none', ...
    'twin_ev',         {[0,0.4]}, ...
    'stat_to_use',     {'ccMeanDiff'}, ...
    'nTrialsMin',      20, ...
    'nTrialsToKeep',   false, ...
    'minTrialsPerCond',10, ...
    'minTrialsPerCombo',2, ...
    'statCV',          '', ...
    'pthresh',         0.05, ...
    'nPseudoSessions', 199 ...
    );

PETH_struct = get_taskPETH_andTuning(datpath,params); 

```

## Data

This repository does not include raw data. 
Analysis scripts expect data to be locally accessible under 'Y:\Subjects'


## Authors

Samuel Picard
