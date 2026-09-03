function [twin_ev, condVals, T] = get_twins(evnt, trialTypeField)
%GET_TWINS Return default response windows and condition values.
%
% Extend as needed.

T = []; % Leave empty unless using a known fixed time base.

switch evnt

    case 'stimOn'

        switch trialTypeField
            case 'contrastDiff'
                twin_ev = [0.2, 0.6];
                condVals = ...
                    [-1,-.5,-.25,-.125,-.0625,0, ...
                     .0625,.125,.25,.5,1];

            case 'probabilityLeft'
                twin_ev = [-0.5, 0];
                condVals = [0.2, 0.8];

            otherwise
                twin_ev = [0, 0.4];
                condVals = [];
        end


    case 'choiceMovement'

        twin_ev = [0, 0.4];

        switch trialTypeField
            case 'choice'
                condVals = [-1, 1];

            otherwise
                condVals = [];
        end


    case 'feedback'

        twin_ev = [0.2, 0.6];

        switch trialTypeField
            case 'feedbackType'
                condVals = [-1, 1];

            otherwise
                condVals = [];
        end


    case {'valveOn', 'toneOn', 'noiseOn'}

        % Passive stimulus response window.
        twin_ev = [0, 0.4];

        % Passive PETHs have no trial-condition dimension.
        condVals = [];


    otherwise

        twin_ev = [0, 0.4];
        condVals = [];

end

end