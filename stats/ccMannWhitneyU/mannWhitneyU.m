 
function numer = mannWhitneyU(x,y)
% function numer = mannWhitneyU(x,y)
%
% numer is the number of instances for which x>y, of all possible
% comparisons. Divide by nx*ny for mannWhitney u statistic.
%
% x and y are vectors
%
% Samuel Picard 2024 (simplified from steinmetz-et-al-2019)

nx = numel(x);
ny = numel(y);

t = tiedrank([x(:); y(:)]);

% if nx==1
%     numer = t(:)';
% else
%     numer = sum(t,1);
% end

R1 = sum(t(1:nx));
numer = R1-nx*(nx+1)/2;

%R2 = sum(t(nx+1:end));
%numer = min([R1-nx*(nx+1)/2,R2-ny*(ny+1)/2]);
