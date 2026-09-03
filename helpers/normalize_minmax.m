function sig_norm = normalize_minmax(sig,minvals,maxvals)

%normalize a n x T signal matrix by n minvals and maxvals (usually
%percentiles)

if size(sig,1)==1
    sig_norm = (sig - minvals) ./ (maxvals-minvals);
else
    if ismatrix(sig)
        sig_norm = (sig - repmat(minvals,[1,size(sig,2)])) ./ repmat((maxvals-minvals),[1,size(sig,2)]);
    elseif ndims(sig)==3
        sig_norm = (sig - repmat(minvals,[1,size(sig,2),size(sig,3)])) ./ repmat((maxvals-minvals),[1,size(sig,2),size(sig,3)]);
    end
end