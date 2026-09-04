function params_all = normalizeParams(params_all, nPseudosDefault)
% Add missing fields and reconcile legacy aliases.
params_all = ensureField(params_all, 'activity_type', 'deconv');
params_all = ensureField(params_all, 'trialTypeVals', []);
params_all = ensureField(params_all, 'trialTypeFilter', '');
params_all = ensureField(params_all, 'twin_bl', 'none');
params_all = ensureField(params_all, 'nTrialsMin', 20);
params_all = ensureField(params_all, 'nTrialsToKeep', false);
params_all = ensureField(params_all, 'minTrialsPerCond', 10);
params_all = ensureField(params_all, 'minTrialsPerCombo', 2);
params_all = ensureField(params_all, 'stat_to_use', 'ccu');
params_all = ensureField(params_all, 'nComp', 1);
params_all = ensureField(params_all, 'pthresh', 0.05);
params_all = ensureField(params_all, 'nPseudoSessions', nPseudosDefault);

if ~isfield(params_all, 'condFields')
    for i = 1:numel(params_all)
        if isfield(params_all, 'cndFields')
            params_all(i).condFields = params_all(i).cndFields;
        else
            params_all(i).condFields = '';
        end
    end
else
    for i = 1:numel(params_all)
        if isempty(params_all(i).condFields)
            params_all(i).condFields = '';
        end
    end
end

if ~isfield(params_all, 'statCV')
    for i = 1:numel(params_all)
        if isfield(params_all, 'cv') && (ischar(params_all(i).cv) || isstring(params_all(i).cv))
            params_all(i).statCV = char(params_all(i).cv);
        else
            params_all(i).statCV = '';
        end
    end
else
    for i = 1:numel(params_all)
        if isempty(params_all(i).statCV), params_all(i).statCV = ''; end
    end
end
end


function S = ensureField(S, name, defaultValue)
if ~isfield(S, name)
    for i = 1:numel(S)
        S(i).(name) = defaultValue;
    end
else
    for i = 1:numel(S)
        if isempty(S(i).(name))
            S(i).(name) = defaultValue;
        end
    end
end
end
