function [saturatedRatio, maxImg, timeSaturatedImg] = ...
    getSaturatedPixelRatio(filePath, H, W, nFramesToLoad)
% getSaturatedPixelRatio
%
% Quantifies saturation in a motion-registered imaging movie.
%
% The movie is assumed to be:
%   - signed 16-bit integer (int16)
%   - little endian
%   - stored as H x W x N
%
% A common saturation cap is assumed across all pixels. A pixel is
% classified as saturated if it reaches the global maximum value in at
% least two frames.
%
% INPUTS
%   filePath        - path to imaging.frames_motionRegistered.bin
%   H               - image height (default: 512)
%   W               - image width  (default: 512)
%   nFramesToLoad   - number of initial frames to analyse (default: 2000)
%
% OUTPUTS
%   saturatedRatio    - fraction of pixels that reach the global maximum
%                       in at least two frames
%
%   maxImg            - H x W image containing the maximum value reached
%                       by each pixel
%
%   timeSaturatedImg  - H x W image containing the fraction of analysed
%                       frames in which each pixel equals the global
%                       saturation value

if nargin < 2 || isempty(H)
    H = 512;
end

if nargin < 3 || isempty(W)
    W = 512;
end

if nargin < 4 || isempty(nFramesToLoad)
    nFramesToLoad = 2000;
end

%% Load movie

fid = fopen(filePath, 'r', 'ieee-le');

if fid == -1
    error('Could not open file: %s', filePath);
end

cleanupObj = onCleanup(@() fclose(fid));

movie = fread(fid, [H*W, nFramesToLoad], '*int16');

nFramesRead = size(movie, 2);

if nFramesRead < nFramesToLoad
    warning('Only %d of %d requested frames available in %s.', ...
        nFramesRead, nFramesToLoad, filePath);
end

if nFramesRead < 2
    saturatedRatio = NaN;
    maxImg = [];
    timeSaturatedImg = [];
    warning('Fewer than 2 frames available in %s.', filePath);
    return
end

movie = reshape(movie, H, W, nFramesRead);

% Maximum value reached by each pixel
maxImg = max(movie, [], 3);

% Common maximum across the entire movie
globalMax = max(maxImg, [], 'all');

% Count frames where each pixel reaches this maximum
nSaturatedFrames = sum(movie == globalMax, 3);

% Count pixels reaching the cap at least twice
saturatedPixels = nSaturatedFrames >= 2;
saturatedRatio = mean(saturatedPixels, 'all');

% Fraction of frames each pixel spends at the cap
timeSaturatedImg = nSaturatedFrames / nFramesToLoad;

%% Report

fprintf('%s\n', filePath);
fprintf('  Global maximum: %d\n', globalMax);
fprintf('  %d / %d pixels saturated (%.3f%%)\n', ...
    sum(saturatedPixels(:)), numel(saturatedPixels), ...
    100*saturatedRatio);

end