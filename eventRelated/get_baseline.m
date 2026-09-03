function bl = get_baseline(sig,frameTimes,trialsTL,evnt,twin_bl)

evnt_bl = evnt;

if nargin<5
    twin_bl = [-0.4,0];
    evnt_bl = "stimOn_times";
end

Fs = 1/median(diff(frameTimes));

frwin_bl = round(Fs*twin_bl);

numcells = size(sig,1);
numtrials = size(trialsTL,1);

%extract raw traces on each trial
bl = nan(numcells,numtrials);
for iTrial = 1:numtrials
    bl_t = trialsTL{iTrial,evnt_bl};
    if ~isnan(bl_t)
        [~,bl_fr] = min(abs(frameTimes-bl_t)); %get 2p frame corresponding to the baseline time
        bl(:,iTrial) = nanmean(sig(:,(frwin_bl(1):frwin_bl(end))+bl_fr),2);
    else
        warning(['Trial nr. ',num2str(iTrial),' has no valid timestamp! Values set to NaN']);
    end
end

end