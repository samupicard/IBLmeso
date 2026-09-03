function acrList = regionTokenMap(token)
% regionTokenMap: expand special tokens into a list of acronyms
% Return {} if token is not a special token.

switch token
    case 'HVAs'
        acrList = {'VISam','VISl','VISm','VISpo','VISpl'};
    case 'PPC'
        acrList = {'VISa','VISrl'};
    otherwise
        acrList = {};
end
end