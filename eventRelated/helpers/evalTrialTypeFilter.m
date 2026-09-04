function mask = evalTrialTypeFilter(trialsT, trialTypeFilter)
% Evaluate trialTypeFilter as a logical expression over columns of trialsT.
%
% Examples:
%   'contrastDiff~=0'
%   'choice==-1'
%   'contrastDiff==0 && feedbackType==-1'
%   'abs(contrastDiff)>=0.25 & choice==-1'
%   'choiceMovement_times - stimOn_times > 0.4'

nTrials = height(trialsT);

if isempty(trialTypeFilter) || strlength(string(trialTypeFilter)) == 0
    mask = true(nTrials,1);
    return
end

expr = char(string(trialTypeFilter));

expr = strrep(expr, '&&', '&');
expr = strrep(expr, '||', '|');

varNames = trialsT.Properties.VariableNames;

for i = 1:numel(varNames)
    v = varNames{i};

    if isvarname(v)
        val = trialsT.(v);

        if isrow(val)
            val = val(:);
        end

        eval(sprintf('%s = val;', v));
    end
end

try
    mask = eval(expr);
catch ME
    error('Could not evaluate trialTypeFilter "%s": %s', trialTypeFilter, ME.message);
end

if ~islogical(mask)
    mask = logical(mask);
end

mask = mask(:);

if numel(mask) ~= nTrials
    error('trialTypeFilter "%s" did not return one logical value per trial.', trialTypeFilter);
end
end