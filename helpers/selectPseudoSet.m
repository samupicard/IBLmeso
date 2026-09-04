function pseudoSessions = selectPseudoSet(S, trialTypeField)
pseudoSessions = [];
if strncmpi(trialTypeField,'contrast',5) || strcmpi(trialTypeField,'stimSide')
    nm = 'pseudoSessions_contrast';
elseif strncmpi(trialTypeField,'choice',5) || strncmpi(trialTypeField,'movement',5)
    nm = 'pseudoSessions_choice';
elseif strncmpi(trialTypeField,'feedback',5)
    nm = 'pseudoSessions_feedback';
elseif strncmpi(trialTypeField,'probability',5) || strncmpi(trialTypeField,'prior',5)
    nm = 'pseudoSessions_bias';
else
    return
end
if isfield(S,nm), pseudoSessions = S.(nm); end
end

