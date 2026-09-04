function cmap = choosePethColormap(pethType)

[~, subtype] = splitPethType(pethType);

switch subtype
    case 'contrastDiff'
        cmap = brewermap(256, '*RdBu');

    case 'choice'
        %cmap = brewermap(256, '*PiYG');
        cmap = brewermap(256, '*RdBu');

    case 'feedbackType'
        %cmap = brewermap(256, 'PuOr');
        cmap = brewermap(256, '*RdBu');

    case 'probabilityLeft'
        %cmap = brewermap(356, '*RdBu'); cmap = cmap(51:end-50,:);
        cmap = brewermap(256, '*RdBu');

    otherwise
        cmap = brewermap(256, '*RdBu'); % fallback
end

end