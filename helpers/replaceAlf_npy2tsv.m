% Root folder containing many subfolders
rootDir = 'Y:\Subjects';

attribute_name = 'taskTunedP';

%Column headers
col_names = {
    'ccu_stimOn_0to400_stimSide'
    'ccu_choiceMovement_-400to0_choice'
    'ccu_feedback_0to400_feedbackType'
    'ccu_stimOn_-500to-100_probabilityLeft'
};

% stat_names = {
%     'signrank_stimOn_vs_baseline'
%     'signrank_choiceMovement_vs_baseline_choice'
%     'ranksum_choiceMovement_choice'
%     'signrank_feedback_vs_baseline_feedbackType'
%     'ranksum_feedback_feedbackType'
%     'ranksum_baseline_block_probabilityLeft'
% };


%Find all .npy files recursively
D = 5;
patternParts = repmat({'*'}, 1, D);
files = dir(fullfile(rootDir, patternParts{:}, ['mpciROIs.',attribute_name,'.npy']));

fprintf('Found %d files.\n', numel(files));

for i = 1:numel(files)
    
    folderPath = files(i).folder;
    npyPath = fullfile(folderPath, files(i).name);
    
    fprintf('Processing: %s\n', npyPath);
    
    % Load numpy array (requires npy-matlab toolbox)
    data = readNPY(npyPath);
    
    % Sanity check
    if size(data,2) ~= length(col_names)
        warning('File %s does not have %d columns. Skipping.', npyPath, length(col_names));
        continue;
    end
    
    % Create output TSV file
    tsvPath = fullfile(folderPath, ['mpciROIs.',attribute_name,'.tsv']);
    
    T = array2table(data, 'VariableNames', col_names);
    writetable(T, tsvPath, 'FileType', 'text', 'Delimiter', '\t');
    
    % Remove npy file
    delete(npyPath);
    
    % ---- Remove names.tsv if it exists ----
    namesFilePath = fullfile(folderPath, [attribute_name,'.names.tsv']);
    if exist(namesFilePath, 'file')
        delete(namesFilePath);
        fprintf('Deleted: %s\n', namesFilePath);
    end
    
end

fprintf('Done.\n');