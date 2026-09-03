function contrastDiff = getContrastDiff(contrastLeft,contrastRight)

%computes an array of contrast differences based on two arrays of contrasts
% (left and right)

trials_contrastRight = contrastRight;
trials_contrastRight(isnan(trials_contrastRight))=0;

trials_contrastLeft = contrastLeft;
trials_contrastLeft(isnan(trials_contrastLeft))=0;

contrastDiff = trials_contrastRight-trials_contrastLeft;