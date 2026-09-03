function cluster_ids = read_uid_csv(path_to_csv)

%helper function for reading a csv containing strings such as
%*.clusterUIDs.csv files, when this is a single column with some empty rows. 
%Outputs an array of strings with the same number of rows as the original csv, 
%leaving empty rows as empty strings.

% Read entire file as lines of text
lines = readlines(path_to_csv);

lines = lines(1:end-1); %last line is blank and should always be removed for some strange reason

% Ensure all lines are strings and preserve empty ones
cluster_ids = strings(size(lines));
for i = 1:numel(lines)
    if strlength(strtrim(lines(i))) == 0
        cluster_ids(i) = "";  % explicitly empty
    else
        cluster_ids(i) = lines(i);
    end
end