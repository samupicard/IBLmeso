function ttl = prettyPethTitle(pethType)
% Return a 2-line title as {topLine; bottomLine}

switch char(string(pethType))
    case 'stimOn_contrastDiff'
        ttl = {'STIM SIDE'; 'stimOn'};
    case 'choiceMovement_choice'
        ttl = {'CHOICE'; 'choiceMov'};
    case 'feedback_feedbackType'
        ttl = {'FEEDBACK'; 'feedback'};
    case 'stimOn_probabilityLeft'
        ttl = {'BLOCK'; 'stimOn'};
    otherwise
        % fallback: subtype on top, event on bottom
        [evnt, subtype] = splitPethType(pethType);

        % optional light prettifying for fallback
        topLine = upper(regexprep(subtype, '([a-z])([A-Z])', '$1 $2'));
        botLine = evnt;

        ttl = {topLine; botLine};
    end
end