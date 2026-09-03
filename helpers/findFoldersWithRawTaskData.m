function matchingParents = findFoldersWithRawTaskData(rootDir,nLevels)
% Efficiently finds folders exactly 4 levels down that contain 'raw_task_data*' subfolders

if nargin<2
    nLevels=4;
end

getNames = @(d) {d.name};
isVisibleDir = @(d) [d.isdir] & ~startsWith(getNames(d), '.');

lvl1 = dir(rootDir);
if strcmp(rootDir,'Y:\Subjects')
    lvl1 = lvl1(isVisibleDir(lvl1) & startsWith(getNames(lvl1), 'SP'));
    nLevels = 4;
elseif contains(rootDir,'SP')
    nLevels = 3;
    lvl1 = lvl1(isVisibleDir(lvl1));
else
    lvl1 = lvl1(isVisibleDir(lvl1));
end
matchingParents = {};

% Step through nLevels down
switch nLevels
    
    case 3
        
        for d1 = lvl1'
            lvl2 = dir(fullfile(rootDir, d1.name)); lvl2 = lvl2(isVisibleDir(lvl2));
            for d2 = lvl2'
                p2 = fullfile(rootDir, d1.name, d2.name);
                lvl3 = dir(fullfile(p2, '*')); lvl3 = lvl3(isVisibleDir(lvl3));
                if any(startsWith({lvl3.name}, 'raw_task_data'))
                    matchingParents{end+1} = p2; %#ok<AGROW>
                end
            end
        end
        
    case 4
        
        for d1 = lvl1'
            lvl2 = dir(fullfile(rootDir, d1.name)); lvl2 = lvl2(isVisibleDir(lvl2));
            for d2 = lvl2'
                lvl3 = dir(fullfile(rootDir, d1.name, d2.name)); lvl3 = lvl3(isVisibleDir(lvl3));
                for d3 = lvl3'
                    p3 = fullfile(rootDir, d1.name, d2.name, d3.name);
                    lvl4 = dir(fullfile(p3, '*')); lvl4 = lvl4(isVisibleDir(lvl4));
                    if any(startsWith({lvl4.name}, 'raw_task_data'))
                        matchingParents{end+1} = p3; %#ok<AGROW>
                    end
                end
            end
        end
        
end