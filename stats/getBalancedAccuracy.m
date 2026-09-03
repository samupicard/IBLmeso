function bacc = getBalancedAccuracy(y, yPred)

%get balanced accuracy score of a binary vector of predictions against
%a binary vector of data. 
%
% Balanced Accuracy = (sensitivity + specificity) / 2
% where
% Sensitivity = True Positives / All Positive predictions
% Specificity = True Negatives / All Negative predictions
%
% Samuel Picard

TP = sum(yPred & yPred==y); %true positives
TN = sum(~yPred & yPred==y); %true negatives
FP = sum(yPred & yPred~=y); %false positives
FN = sum(~yPred & yPred~=y); %false negatives

Sensitivity = TP / (TP+FN);
Specificity = TN / (TN+FP);

bacc =(Sensitivity + Specificity) / 2;
        