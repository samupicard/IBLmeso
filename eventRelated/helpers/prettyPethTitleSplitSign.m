function ttl = prettyPethTitleSplitSign(pethType, signChar)
% prettyPethTitleSplitSign
%
% Return a 2-line title for splitSign mode, where signChar is '+' or '-'.
%
% Examples:
%   stimOn_contrastDiff (+) -> {'STIM L'; 'stimOn'}
%   stimOn_contrastDiff (-) -> {'STIM R'; 'stimOn'}

pethType = char(string(pethType));
signChar = char(string(signChar));

switch pethType
    case 'stimOn_contrastDiff'
        switch signChar
            case '+'
                ttl = {'stim L'; 'stimOn'};
            case '-'
                ttl = {'stim R'; 'stimOn'};
            otherwise
                error('signChar must be ''+'' or ''-''.');
        end
    case 'choiceMovement_choice'
        switch signChar
            case '+'
                ttl = {'choice L'; 'choiceMov'};
            case '-'
                ttl = {'choice R'; 'choiceMov'};
            otherwise
                error('signChar must be ''+'' or ''-''.');
        end

    case 'feedback_feedbackType'
        switch signChar
            case '+'
                ttl = {'fb +'; 'feedback'};
            case '-'
                ttl = {'fb -'; 'feedback'};
            otherwise
                error('signChar must be ''+'' or ''-''.');
        end

    case 'stimOn_probabilityLeft'
        switch signChar
            case '+'
                ttl = {'block L'; 'stimOn'};
            case '-'
                ttl = {'block R'; 'stimOn'};
            otherwise
                error('signChar must be ''+'' or ''-''.');
        end

    otherwise
        % fallback
        base = prettyPethTitle(pethType);
        ttl = {sprintf('%s %s', base{1}, signChar); base{2}};
    end
end