function [Evk,Dat,T] = get_evoked_full(sig,frameTimes,timeShifts,trialsTL,evnt,twin_ev,twin_bl,twin_all)

if nargin<7; 
    twin_all = [-1.6,3]; 
end


frame_intervals = diff(frameTimes);
Fs = 1/median(frame_intervals); %original sampling frequency in Hz
if Fs<20
    Fs_us = 20; %upsampling frequency in Hz
end

jump_idxs = [find(frame_intervals>2*Fs) length(frameTimes)];
frameTimes_us = [];
block_start = 1;
for iBlock = 1:length(jump_idxs)
    block_end = jump_idxs(iBlock);
    frameTimes_us = [frameTimes_us frameTimes(block_start):1/Fs_us:frameTimes(block_end)];
    block_start = block_end+1;
end

sig_us = nan(size(sig,1),length(frameTimes_us));
frameTimes_us_shifts = nan(size(sig,1),length(frameTimes_us));
tic
for i = 1:size(sig,1)
    sig_us(i,:) = interp1(frameTimes,sig(i,:),frameTimes_us);
    frameTimes_us_shifts(i,:) = frameTimes_us+timeShifts(i);
end
toc

frwin_all = round(Fs_us*twin_all);
frwin_bl = round(Fs_us*twin_bl);
frwin_ev = round(Fs_us*twin_ev);

T = [frwin_all(1):frwin_all(end)]/Fs_us;

numcells = size(sig_us,1);
numframes = size(T,2);
numtrials = size(trialsTL,1);

%extract raw traces on each trial
Dat = nan(numcells,numtrials,numframes);
Evk = nan(numcells,numtrials);
for iTrial = 1:numtrials
    evnt_t = trialsTL{iTrial,evnt}; %look up event timestamp
    if ~isnan(evnt_t)
        tic; 
        for i=1:size(sig,1)
        [~,evnt_fr] = min(abs(frameTimes_us+timeShifts(i)-evnt_t)); %get 2p frame corresponding to the event time
        bl = nanmean(sig_us(i,(frwin_bl(1):frwin_bl(end))+evnt_fr),2);
        Dat(i,iTrial,:) = sig_us(i,(frwin_all(1):frwin_all(end))+evnt_fr);
        %Dat(:,iOri,iRep,:) = sigS(:,[frwin_all(1):frwin_all(end)]+evnt_frame); %smoothed
        Evk(i,iTrial) = nanmean(sig_us(i,(frwin_ev(1):frwin_ev(end))+evnt_fr),2)-bl;
        %Evk(:,iTrial) = nanmean(sig_us(:,[frwin_ev(1):frwin_ev(end)]+evnt_frame),2); %without baseline correction
        end
        toc;
    else
        warning(['Trial nr. ',num2str(iTrial),' has no valid timestamp! Values set to NaN']);
    end
end

end