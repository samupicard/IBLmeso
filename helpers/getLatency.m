function trialsT_lats = getLatency(trialsT, evnt)
%getLatency Subtracts reference column values from all columns ending in 'times'.
%
%   trialsT_lats = getLatency(T, refColName) subtracts the values in the
%   reference column (refColName) from each column in trialsT whose name ends
%   with 'times'. The result is returned in trialsT_lats.
%
%   Inputs:
%     trialsT   - Input trials table (nTrials x nCols).
%     evnt      - Name (string or char) of the reference column to subtract.
%
%   Output:
%     trialsT_lats  - Output table with updated time fields.
%
% Samuel Picard

s = '_times'; 
s_new = '_lats';

%append event name with _times
N = min(length(evnt),length(s));
if strcmp(evnt(end-(N-1):end), s)
    evnt_nm = evnt;
    evnt = evnt_nm(1:end-length(s));
else
    evnt_nm = [evnt s];
end

%get all column names
colNames = trialsT.Properties.VariableNames;
isTimeCol = endsWith(colNames,s);

%make new table copied from the first
trialsT_lats = trialsT;

%timeCols = colNames(endsWith(colNames, s));
%trialsT_lats = table();

for i = 1:numel(colNames)
    %col = timeCols{i};
    if isTimeCol(i)
        col = colNames{i};
        col_new = [col(1:end-N) s_new];
        trialsT_lats.Properties.VariableNames{i} = col_new;
        trialsT_lats.(col_new) = trialsT.(col) - trialsT.(evnt_nm);
    end
end

end
