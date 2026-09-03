function logL = getLogL(y, p)

%get log likelihood of a predicted probability of a binary event, against
%a binary vector of data (Bernoulli log-likelihood). 
%
% for each event t,
% LogL = sum( y(t)*log(p(t)) + (1-y(t))*log(1-p(t)) )
%
% Samuel Picard

%make sure y and p have the correct sizes
if size(y,1)>=1 && size(y,2)==1
    y = y';
end
if size(p,2)>=1 && size(p,1)==1
    p = p';
end
%size(p)
%make sure we have a binary vector
classes = unique(y);
if length(classes)>2
    error('There are more than 2 classes in y.')
end
if ~all(ismember(y,[0,1]))
    y = double(y==classes(end));
end

%make p ranges between 0 and 1
if any(p)<0 || any(p)>1
    error('Probabilities p are not between 0 and 1.')
end

logL = y*log(p) + (1-y)*log(1-p);
        