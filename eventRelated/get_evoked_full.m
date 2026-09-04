function [Evk,Dat,T,badTrials,badTrialsEvk] = get_evoked_full(sig,frameTimes,badframes,trialsTL,evnt,twin_ev,twin_bl,twin_all,timeShifts)
%GET_EVOKED_FULL Return time-averaged responses and full event-aligned data.
%
% [Evk,Dat,T,badTrials,badTrialsEvk] = get_evoked_full(...)
%
% Evk           : nCells x nTrials mean response in twin_ev (baseline
%                 subtracted when twin_bl is numeric).
% Dat           : nCells x nTrials x nTimepoints full PETH in twin_all.
% T             : relative time vector corresponding to Dat.
% badTrials     : trial invalid for the FULL PETH (twin_all). 
% badTrialsEvk  : trial invalid for the event-response calculation only
%                 (twin_ev and, when used, twin_bl). This is useful when
%                 Evk is used for tuning statistics: a trial can have a
%                 valid Evk even when twin_all extends outside valid data.
%
% timeShifts can be empty, one scalar shift per cell, or a matrix matching
% sig (cell x frame).

if nargin < 9
    timeShifts = [];
end
if nargin < 8
    twin_all = [-0.9,3.3];
end

if ischar(twin_bl) || (isstring(twin_bl) && isscalar(twin_bl))
    if strcmp(char(twin_bl),'none')
        bl_flag = false;
    else
        warning('Unclear whether to baseline correct. Not doing baseline correction!');
        bl_flag = false;
    end
else
    bl_flag = true;
end

Fs = 1/median(diff(frameTimes));

% temporal resolution guardrail
dt = median(diff(frameTimes));
evWin = abs(diff(twin_ev));
minSamplesEv = evWin / dt;
if minSamplesEv < 1 - 1e-5
    warning(['Event window twin_ev = [%.4f %.4f] is only %.4fs wide, ' ...
             'but imaging dt is %.4fs (Fs = %.2f Hz). ' ...
             'This gives only %.2f expected samples per window. ' ...
             'Consider resampling before calling get_evoked_full, or using a wider twin_ev.'], ...
             twin_ev(1), twin_ev(2), evWin, dt, Fs, minSamplesEv);
end

if bl_flag
    blWin = abs(diff(twin_bl));
    minSamplesBl = blWin / dt;
    if minSamplesBl < 1 - 1e-5
        warning(['Baseline window twin_bl = [%.4f %.4f] is only %.4fs wide, ' ...
                 'but imaging dt is %.4fs (Fs = %.2f Hz). ' ...
                 'This gives only %.2f expected samples per window. ' ...
                 'Consider resampling before calling get_evoked_full, or using a wider twin_bl.'], ...
                 twin_bl(1), twin_bl(2), blWin, dt, Fs, minSamplesBl);
    end
end

numcells = size(sig,1);
numtrials = size(trialsTL,1);
nSigFrames = size(sig,2);

frameTimes = frameTimes(:)';

frwin_all = round(Fs*twin_all);
T = (frwin_all(1):frwin_all(end)) / Fs;
numframes = numel(T);

Evk = nan(numcells,numtrials);
Dat = nan(numcells,numtrials,numframes);
badTrials = false(1,numtrials);
badTrialsEvk = false(1,numtrials);

badMask = false(1,nSigFrames);
badframes = badframes(badframes >= 1 & badframes <= nSigFrames);
badMask(badframes) = true;

useTimeShifts = ~isempty(timeShifts);

if useTimeShifts
    if isvector(timeShifts)
        timeShifts = timeShifts(:);
        if numel(timeShifts) ~= numcells
            error('get_evoked_full:InvalidTimeShifts', ...
                'Vector timeShifts must contain one value per cell (%d).', numcells);
        end
    elseif ~isequal(size(timeShifts), size(sig))
        error('get_evoked_full:InvalidTimeShifts', ...
            'Matrix timeShifts must have the same size as sig.');
    end

    minShift = min(timeShifts(:));
    maxShift = max(timeShifts(:));
end

nCharsPrinted = 0;
fprintf('trial ');

for iTrial = 1:numtrials

    fprintf(repmat('\b',1,nCharsPrinted))
    nCharsPrinted = fprintf('%d/%d ..',iTrial,numtrials);

    try
        evnt_t = getEventTimestamp(trialsTL, iTrial, evnt);

        if isempty(evnt_t) || ~isscalar(evnt_t) || ~isfinite(evnt_t) || ...
                min(abs(frameTimes - evnt_t)) > (1/Fs)
            badTrialsEvk(iTrial) = true;
            badTrials(iTrial) = true;
            continue
        end

        if ~useTimeShifts
            %% Event response validity and value
            evIdx = find(frameTimes >= evnt_t + twin_ev(1) & ...
                         frameTimes <= evnt_t + twin_ev(2));

            if isempty(evIdx) || any(badMask(evIdx))
                badTrialsEvk(iTrial) = true;
                badTrials(iTrial) = true;
                continue
            end

            ev = mean(sig(:,evIdx),2,'omitnan');

            if bl_flag
                blIdx = find(frameTimes >= evnt_t + twin_bl(1) & ...
                             frameTimes <= evnt_t + twin_bl(2));

                if isempty(blIdx) || any(badMask(blIdx))
                    badTrialsEvk(iTrial) = true;
                    badTrials(iTrial) = true;
                    continue
                end

                bl = mean(sig(:,blIdx),2,'omitnan');
                Evk(:,iTrial) = ev - bl;
            else
                bl = [];
                Evk(:,iTrial) = ev;
            end

            %% Full PETH validity and value
            [~,evnt_fr] = min(abs(frameTimes - evnt_t));
            idx_all = (frwin_all(1):frwin_all(end)) + evnt_fr;

            if any(idx_all < 1) || any(idx_all > nSigFrames) || any(badMask(idx_all))
                badTrials(iTrial) = true;
                continue
            end

            if bl_flag
                Dat(:,iTrial,:) = sig(:,idx_all) - bl;
            else
                Dat(:,iTrial,:) = sig(:,idx_all);
            end

        else
            %% Event response validity and value, allowing cell-specific shifts
            evCand = find(frameTimes >= evnt_t + twin_ev(1) - maxShift & ...
                          frameTimes <= evnt_t + twin_ev(2) - minShift);

            if isempty(evCand)
                badTrialsEvk(iTrial) = true;
                badTrials(iTrial) = true;
                continue
            end

            if isvector(timeShifts)
                evTimes = frameTimes(evCand) + timeShifts;
            else
                evTimes = frameTimes(evCand) + timeShifts(:,evCand);
            end

            evMask = evTimes >= evnt_t + twin_ev(1) & ...
                     evTimes <= evnt_t + twin_ev(2);
            evBad = evMask & badMask(evCand);

            evVals = sig(:,evCand);
            validEvSamples = evMask & ~isnan(evVals);
            evVals(~validEvSamples) = 0;
            evCount = sum(validEvSamples,2);
            evSum = sum(evVals,2);

            badCellEvk = evCount == 0 | any(evBad,2);
            ev = evSum ./ evCount;

            if bl_flag
                blCand = find(frameTimes >= evnt_t + twin_bl(1) - maxShift & ...
                              frameTimes <= evnt_t + twin_bl(2) - minShift);

                if isempty(blCand)
                    badTrialsEvk(iTrial) = true;
                    badTrials(iTrial) = true;
                    continue
                end

                if isvector(timeShifts)
                    blTimes = frameTimes(blCand) + timeShifts;
                else
                    blTimes = frameTimes(blCand) + timeShifts(:,blCand);
                end

                blMask = blTimes >= evnt_t + twin_bl(1) & ...
                         blTimes <= evnt_t + twin_bl(2);
                blBad = blMask & badMask(blCand);

                blVals = sig(:,blCand);
                validBlSamples = blMask & ~isnan(blVals);
                blVals(~validBlSamples) = 0;
                blCount = sum(validBlSamples,2);
                blSum = sum(blVals,2);

                badCellEvk = badCellEvk | blCount == 0 | any(blBad,2);
                bl = blSum ./ blCount;

                goodCellEvk = ~badCellEvk;
                Evk(goodCellEvk,iTrial) = ev(goodCellEvk) - bl(goodCellEvk);
            else
                bl = [];
                goodCellEvk = ~badCellEvk;
                Evk(goodCellEvk,iTrial) = ev(goodCellEvk);
            end

            if any(badCellEvk)
                badTrialsEvk(iTrial) = true;
            end

            %% Full PETH validity and value
            % Find the event-centre frame for each cell without constructing
            % an nCells x nFrames matrix. This was previously the dominant
            % cost when timeShifts was supplied.
            if isvector(timeShifts)
                % One constant time shift per cell. The required unshifted
                % frame time is simply evnt_t - timeShift. interp1 performs
                % the nearest-frame lookup in 1-D.
                targetTimes = evnt_t - timeShifts;
                evnt_fr = interp1(frameTimes, 1:nSigFrames, targetTimes, ...
                    'nearest', NaN);
            else
                % Frame-dependent shifts cannot be inverted analytically.
                % Restrict the exhaustive search to frames that could
                % possibly align to this event given the global shift range.
                centreCand = find(frameTimes >= evnt_t - maxShift - dt & ...
                                  frameTimes <= evnt_t - minShift + dt);

                if isempty(centreCand)
                    badTrials(iTrial) = true;
                    continue
                end

                shiftedCentreTimes = frameTimes(centreCand) + timeShifts(:,centreCand);
                [~,localIdx] = min(abs(shiftedCentreTimes - evnt_t),[],2);
                evnt_fr = centreCand(localIdx);
            end

            win_all = frwin_all(1):frwin_all(end);
            idx_all = evnt_fr + win_all;

            badCellDat = isnan(evnt_fr) | ...
                         any(idx_all < 1 | idx_all > nSigFrames,2);

            % Only index rows whose whole window is in range. Use direct
            % column-major linear indexing instead of repmat + sub2ind.
            inRangeCells = ~badCellDat;
            if any(inRangeCells)
                rowsInRange = find(inRangeCells);
                idxGood = idx_all(inRangeCells,:);
                lin_all = rowsInRange + (idxGood - 1) * numcells;
                dat_i = sig(lin_all);

                badCellDat(inRangeCells) = any(badMask(idxGood),2);

                goodCellDat = inRangeCells & ~badCellDat & ~badCellEvk;
                if any(goodCellDat)
                    goodRows = find(goodCellDat);
                    [~,loc] = ismember(goodRows, rowsInRange);
                    if bl_flag
                        Dat(goodRows,iTrial,:) = dat_i(loc,:) - bl(goodRows);
                    else
                        Dat(goodRows,iTrial,:) = dat_i(loc,:);
                    end
                end
            end

            % A baseline/event failure also makes the full PETH unusable for
            % that cell because baseline-corrected Dat cannot be constructed.
            badCellDat = badCellDat | badCellEvk;
            if any(badCellDat)
                badTrials(iTrial) = true;
            end
        end

    catch
        badTrialsEvk(iTrial) = true;
        badTrials(iTrial) = true;
    end
end

fprintf('. Done!\n');

nBadDat = sum(badTrials);
if nBadDat > 0
    warning('%d trials had an invalid full PETH window; Dat values remain NaN for those trials/cells.', nBadDat);
end

nBadEvk = sum(badTrialsEvk);
if nBadEvk > 0
    warning('%d trials had no valid event-response window; Evk values remain NaN for those trials/cells.', nBadEvk);
end

end


function evnt_t = getEventTimestamp(eventData, iTrial, evnt)
%GETEVENTTIMESTAMP Retrieve one event timestamp.
%
% Supports table/timetable variable names or indices, numeric matrices and
% cell arrays.

if istable(eventData) || istimetable(eventData)

    if ischar(evnt) || (isstring(evnt) && isscalar(evnt))
        variableName = char(evnt);

        if ~ismember(variableName,eventData.Properties.VariableNames)
            error('get_evoked_full:UnknownEvent', ...
                'Event variable "%s" is not present in the table.', variableName);
        end

        evnt_t = eventData{iTrial,variableName};

    elseif isnumeric(evnt) && isscalar(evnt)
        evnt_t = eventData{iTrial,evnt};
    else
        error('get_evoked_full:InvalidEventSelector', ...
            'For table inputs, evnt must be a variable name or scalar numeric variable index.');
    end

elseif iscell(eventData)
    evnt_t = eventData{iTrial,evnt};
elseif isnumeric(eventData) || islogical(eventData)
    evnt_t = eventData(iTrial,evnt);
else
    error('get_evoked_full:UnsupportedEventData', ...
        'Unsupported event-data type: %s.',class(eventData));
end

if iscell(evnt_t) && isscalar(evnt_t)
    evnt_t = evnt_t{1};
end
if ischar(evnt_t) || isstring(evnt_t)
    evnt_t = str2double(string(evnt_t));
end
if isnumeric(evnt_t) || islogical(evnt_t)
    evnt_t = double(evnt_t);
end
end
