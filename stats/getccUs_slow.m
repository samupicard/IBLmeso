function Us = getccUs_slow(resps,trialsT,iVariable,trialTypeField)

%get the conditions we want to match
condFields = getBalancedConds(trialTypeField);
conds_to_match = condFields(~strcmp(condFields,trialTypeField));

%find unique conditions
trialsT_conds = trialsT(:,conds_to_match);
[uniqueConds,~,trialCondition] = unique(trialsT_conds,'rows');
    
% %create a shuffle matrix
% numShuffles = 1000;
% [~,shufLabels] = sort(rand(size(trialsT,1),numShuffles,1));
% shufLabels = [(1:size(trialsT,1))' shufLabels];
% 
% %compute combined condition Mann-Whitney U statistic
% [Us, Ups, ccUSummary] = ccUShuf(validResps, iVariable, trialCondition', shufLabels);

Us = nan(1,size(resps,1));
for i=1:size(resps,1)
    [Us(i), ~] = ccU(resps(i,:), iVariable, trialCondition'); %conditions combined Mann-Whitney U-test
end