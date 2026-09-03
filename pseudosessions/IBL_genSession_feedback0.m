function pseudoS = IBL_genSession_feedback0(trialsT)

% IBL_genSession_feedback0 generates new trial feedback for 0-contrast trials, according to biasChoiceWorld heuristics
% we are randomizing the 'correct side' for 0-contrast trials, using draws with probabilityLeft 
% we are keeping everything else the same (block, contrast, choice)
%
% choice=1 means LEFT choice (i.e. the correct choice for LEFT stims)
% choice=-1 means RIGHT choice (i.e. the correct choice for RIGHT stims)

condFields_names = {'probabilityLeft','contrastDiff', 'choice', 'feedbackType'};

pseudoSessT = trialsT(:,condFields_names);

idx = pseudoSessT.contrastDiff == 0;

leftCorrect = rand(sum(idx),1) < pseudoSessT.probabilityLeft(idx); %Bernoulli draws

%reward condition: feedback = +1 when choice == (2*leftCorrect - 1)
correctSide = 2*leftCorrect - 1; 
pseudoSessT.feedbackType(idx) = 2*(pseudoSessT.choice(idx) == correctSide) - 1;

%OLD VERSION
%nTrials = size(trialsT,1);
% for iTrial = 1:nTrials
%     if pseudoSessT{iTrial,'contrastDiff'}==0
%         leftCorrect = rand < pseudoSessT{iTrial,'probabilityLeft'};
%         if leftCorrect && pseudoSessT{iTrial,'choice'}==1 || ~leftCorrect && pseudoSessT{iTrial,'choice'}==-1
%             pseudoSessT.feedbackType(iTrial) = 1;
%         else
%             pseudoSessT.feedbackType(iTrial) = -1;
%         end
%     end
% end

pseudoS = table2struct(pseudoSessT,'ToScalar',true);

end