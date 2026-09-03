function stat_names = infer_stat_names(params_all)
% infer_stat_names
% Build the canonical stat_names strings from params_all.
%
% Convention:
%   <stat>_<event>_<t0ms>to<t1ms>_<trialTypeField-or-special>
%
% Special-case used in your legacy names:
%   contrastDiff comparisons become stimSide or stimSide100 depending on nComp.

    n = numel(params_all);
    stat_names = cell(1, n);

    for i = 1:n
        p = params_all(i);

        % stat
        if isfield(p, 'stat_to_use') && ~isempty(p.stat_to_use)
            stat = char(p.stat_to_use);
        else
            stat = 'stat';
        end

        % event
        evnt = char(p.evnt);

        % time window (ms)
        twin = p.twin_ev;
        t0 = round(1000 * twin(1));
        t1 = round(1000 * twin(2));
        twinStr = sprintf('%dto%d', t0, t1);

        % trialType: default to trialTypeField
        trialType = char(p.trialTypeField);
        
        % special-case mapping for contrastDiff
        if startsWith(trialType, 'contrast', 'IgnoreCase', true)
            % If nComp == 1 -> stimSide100; otherwise stimSide
            if isfield(p,'nComp') && isscalar(p.nComp) && p.nComp == 1
                trialType = 'stimSide100';
            else
                trialType = 'stimSide';
            end
        end

        % If trialTypeFilter isn't empty, append it
        if isfield(p,'trialTypeFilter') && ~isempty(p.trialTypeFilter)
            trialType = [trialType,'_',p.trialTypeFilter];
        end

        % include nTrialsToKeep if it's used
        if p.nTrialsToKeep,
            nTrials = sprintf('_first%d',p.nTrialsToKeep);
        else
            nTrials = '';
        end

        %include cv partition if used
        if isfield(p,'cv') && ~isempty(p.cv)
            cvType = ['_' p.cv];
        else
            cvType = '';
        end

        stat_names{i} = sprintf('%s_%s_%s_%s%s%s', stat, evnt, twinStr, trialType, nTrials, cvType);
    end
end