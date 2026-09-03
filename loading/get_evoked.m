function [Evk,badTrials] = get_evoked( ...
    sig,frameTimes,badframes,trialsTL,evnt,twin_ev,twin_bl,timeShifts)

if nargin < 8
    timeShifts = [];
end

if strcmp(class(twin_bl),'char')
    if strcmp(twin_bl,'none')
        bl_flag = false;
    else
        warning(['Unclear whether to baseline correct. ' ...
                 'Not doing baseline correction!']);
        bl_flag = false;
    end
else
    bl_flag = true;
end

Fs = 1/median(diff(frameTimes));

% Temporal resolution guardrail
dt = median(diff(frameTimes));
evWin = abs(diff(twin_ev));
minSamplesEv = evWin / dt;

if minSamplesEv < 1 - 1e-5
    warning(['Event window twin_ev = [%.4f %.4f] is only %.4fs wide, ' ...
             'but imaging dt is %.4fs (Fs = %.2f Hz). ' ...
             'This gives only %.2f expected samples per window. ' ...
             'Consider resampling before calling get_evoked, ' ...
             'or using a wider twin_ev.'], ...
             twin_ev(1), twin_ev(2), evWin, dt, Fs, minSamplesEv);
end

if bl_flag
    blWin = abs(diff(twin_bl));
    minSamplesBl = blWin / dt;

    if minSamplesBl < 1 - 1e-5
        warning(['Baseline window twin_bl = [%.4f %.4f] is only %.4fs wide, ' ...
                 'but imaging dt is %.4fs (Fs = %.2f Hz). ' ...
                 'This gives only %.2f expected samples per window. ' ...
                 'Consider resampling before calling get_evoked, ' ...
                 'or using a wider twin_bl.'], ...
                 twin_bl(1), twin_bl(2), blWin, dt, Fs, minSamplesBl);
    end
end

numcells = size(sig,1);
numtrials = size(trialsTL,1);
nSigFrames = size(sig,2);

frameTimes = frameTimes(:)';

Evk = nan(numcells,numtrials);
badTrials = false(1,numtrials);

badMask = false(1,nSigFrames);
badMask(badframes) = true;

useTimeShifts = ~isempty(timeShifts);

if useTimeShifts
    if isvector(timeShifts)
        timeShifts = timeShifts(:);
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
        % Supports trials tables, passive Gabor tables, and passive
        % stimulus tables. Row names are ignored.
        evnt_t = getEventTimestamp(trialsTL,iTrial,evnt);

        if isempty(evnt_t) || ...
                ~isscalar(evnt_t) || ...
                ~isfinite(evnt_t) || ...
                min(abs(frameTimes-evnt_t)) > (1/Fs)

            badTrials(iTrial) = true;
            continue
        end

        if ~useTimeShifts

            evIdx = find( ...
                frameTimes >= evnt_t + twin_ev(1) & ...
                frameTimes <= evnt_t + twin_ev(2));

            if isempty(evIdx) || any(badMask(evIdx))
                badTrials(iTrial) = true;
                continue
            end

            ev = mean(sig(:,evIdx),2,'omitnan');

            if bl_flag
                blIdx = find( ...
                    frameTimes >= evnt_t + twin_bl(1) & ...
                    frameTimes <= evnt_t + twin_bl(2));

                if isempty(blIdx) || any(badMask(blIdx))
                    badTrials(iTrial) = true;
                    continue
                end

                bl = mean(sig(:,blIdx),2,'omitnan');
                Evk(:,iTrial) = ev - bl;
            else
                Evk(:,iTrial) = ev;
            end

        else

            % Candidate frames only: shifted frameTimes must fall
            % inside twin_ev.
            evCand = find( ...
                frameTimes >= evnt_t + twin_ev(1) - maxShift & ...
                frameTimes <= evnt_t + twin_ev(2) - minShift);

            if isempty(evCand)
                badTrials(iTrial) = true;
                continue
            end

            if isvector(timeShifts)
                evTimes = frameTimes(evCand) + timeShifts;
            else
                evTimes = frameTimes(evCand) + timeShifts(:,evCand);
            end

            evMask = ...
                evTimes >= evnt_t + twin_ev(1) & ...
                evTimes <= evnt_t + twin_ev(2);

            evBad = evMask & badMask(evCand);

            evVals = sig(:,evCand);
            evVals(~evMask | isnan(evVals)) = 0;

            evCount = sum( ...
                evMask & ~isnan(sig(:,evCand)),2);

            evSum = sum(evVals,2);

            badCell = ...
                evCount == 0 | ...
                any(evBad,2);

            ev = evSum ./ evCount;

            if bl_flag

                blCand = find( ...
                    frameTimes >= evnt_t + twin_bl(1) - maxShift & ...
                    frameTimes <= evnt_t + twin_bl(2) - minShift);

                if isempty(blCand)
                    badTrials(iTrial) = true;
                    continue
                end

                if isvector(timeShifts)
                    blTimes = frameTimes(blCand) + timeShifts;
                else
                    blTimes = ...
                        frameTimes(blCand) + timeShifts(:,blCand);
                end

                blMask = ...
                    blTimes >= evnt_t + twin_bl(1) & ...
                    blTimes <= evnt_t + twin_bl(2);

                blBad = blMask & badMask(blCand);

                blVals = sig(:,blCand);
                blVals(~blMask | isnan(blVals)) = 0;

                blCount = sum( ...
                    blMask & ~isnan(sig(:,blCand)),2);

                blSum = sum(blVals,2);

                badCell = ...
                    badCell | ...
                    blCount == 0 | ...
                    any(blBad,2);

                bl = blSum ./ blCount;

                goodCell = ~badCell;
                Evk(goodCell,iTrial) = ...
                    ev(goodCell) - bl(goodCell);

            else
                goodCell = ~badCell;
                Evk(goodCell,iTrial) = ev(goodCell);
            end

            if any(badCell)
                badTrials(iTrial) = true;
            end
        end

    catch
        badTrials(iTrial) = true;
    end
end

fprintf('. Done!\n');

nantrial_cnt = sum(badTrials);

if nantrial_cnt > 0
    warning([num2str(nantrial_cnt), ...
        ' trials had no valid timestamp or matching imaging frames! ' ...
        'Values set to NaN.']);
end

end


function evnt_t = getEventTimestamp(eventData,iTrial,evnt)
% getEventTimestamp Retrieves one event time from a table or numeric array.
%
% For tables, evnt can be:
%   'stimOn_times'
%   "stimOn_times"
%   'valveOn'
%   'toneOn'
%   'noiseOn'
%   a numeric variable index
%
% Table row names are ignored.

if istable(eventData) || istimetable(eventData)

    if ischar(evnt) || (isstring(evnt) && isscalar(evnt))
        variableName = char(evnt);

        if ~ismember( ...
                variableName, ...
                eventData.Properties.VariableNames)

            error('get_evoked:UnknownEvent', ...
                'Event variable "%s" is not present in the table.', ...
                variableName);
        end

        evnt_t = eventData{iTrial,variableName};

    elseif isnumeric(evnt) && isscalar(evnt)
        evnt_t = eventData{iTrial,evnt};

    else
        error('get_evoked:InvalidEventSelector', ...
            ['For a table input, evnt must be a variable name ' ...
             'or a scalar numeric variable index.']);
    end

else
    % Preserve support for numeric or cell-array inputs.
    evnt_t = eventData{iTrial,evnt};
end

% readtable can occasionally produce a one-element cell.
if iscell(evnt_t) && isscalar(evnt_t)
    evnt_t = evnt_t{1};
end

% Convert scalar strings or character values to numeric timestamps.
if isstring(evnt_t) || ischar(evnt_t)
    evnt_t = str2double(evnt_t);
end

% Ensure numeric types such as single are returned consistently.
if isnumeric(evnt_t) || islogical(evnt_t)
    evnt_t = double(evnt_t);
end

end