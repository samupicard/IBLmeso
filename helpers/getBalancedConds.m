function condFields = getBalancedConds(trialTypeField)

%getBalancedConds defines the marginal conditions to balance out 
% when testing a particular task variable

    if strcmp(trialTypeField,'choice')
        condFields = {'contrastDiff'};
        %condFields = {'probabilityLeft','contrastDiff'};
    elseif strcmp(trialTypeField,'contrastDiff')
        condFields = {'choice'};
        %condFields = {'probabilityLeft','choice'};
    elseif strcmp(trialTypeField,'feedbackType')
        condFields = {'choice'};
        %condFields = {'probabilityLeft','choice'};
    elseif strcmp(trialTypeField,'probabilityLeft')
        condFields = {'contrastDiff','choice'};
    elseif strcmp(trialTypeField,'movement')
        condFields = {'contrastDiff','choice'};
    end

end