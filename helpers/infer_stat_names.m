function stat_names = infer_stat_names(params_all)
%INFER_STAT_NAMES Build canonical stat_names strings from params_all.
%
% Convention:
%   <stat>_<event>_<t0ms>to<t1ms>_<trialTypeField>
%
% Examples:
%   ccu_stimOn_0to400_contrastDiff
%   ccu_stimOn_0to400_stimSide
%   ccu_choiceMovement_-200to200_choice
%   ccu_feedback_0to400_feedbackType
%   ccu_stimOn_-500to-100_probabilityLeft

n = numel(params_all);
stat_names = cell(1,n);

for i = 1:n
    p = params_all(i);

    % statistic
    if isfield(p,'stat_to_use') && ~isempty(p.stat_to_use)
        stat = char(p.stat_to_use);
    else
        stat = 'stat';
    end

    % event
    evnt = char(p.evnt);

    % event-response window in ms
    twin = p.twin_ev;
    t0 = round(1000 * twin(1));
    t1 = round(1000 * twin(2));
    twinStr = sprintf('%dto%d',t0,t1);

    % analysis variable
    trialType = char(p.trialTypeField);

    % optional trial-count restriction
    if isfield(p,'nTrialsToKeep') && p.nTrialsToKeep
        nTrials = sprintf('_first%d',p.nTrialsToKeep);
    else
        nTrials = '';
    end

    % optional CV partition
    if isfield(p,'statCV') && ~isempty(p.statCV)
        cvType = ['_' char(p.statCV)];
    elseif isfield(p,'cv') && ~isempty(p.cv)
        % backwards compatibility
        cvType = ['_' char(p.cv)];
    else
        cvType = '';
    end

    stat_names{i} = sprintf('%s_%s_%s_%s%s%s', ...
        stat, evnt, twinStr, trialType, nTrials, cvType);
end

end

