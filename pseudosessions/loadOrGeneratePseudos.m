function pseudoSets = loadOrGeneratePseudos(datpath, trialsT, nNeeded)
pseudoPath = fullfile(datpath, 'alf', 'pseudoSessions.mat');
regenerate = true;
pseudoSets = struct();

if exist(pseudoPath, 'file')
    S = load(pseudoPath);
    if isfield(S,'pseudoSessions_contrast') && ~isempty(S.pseudoSessions_contrast)
        sameN = numel(S.pseudoSessions_contrast(1).probabilityLeft) == height(trialsT);
        enough = numel(S.pseudoSessions_contrast) >= nNeeded;
        if sameN && enough
            pseudoSets = S;
            regenerate = false;
            fprintf('Loaded existing pseudosessions.\n');
        end
    end
end

if ~regenerate, return; end
fprintf('Generating %d pseudosessions... ', nNeeded);

protoIsBiased = numel(unique(trialsT.probabilityLeft)) == 3;
template = struct('probabilityLeft',[], 'contrastDiff',[], 'choice',[], 'feedbackType',[]);
pseudoSessions_contrast = repmat(template,1,nNeeded);
pseudoSessions_choice   = repmat(template,1,nNeeded);
pseudoSessions_feedback = repmat(template,1,nNeeded);
if protoIsBiased
    pseudoSessions_bias = repmat(template,1,nNeeded);
end

for iP = 1:nNeeded
    if protoIsBiased
        pseudoSessions_bias(iP) = IBL_genSession(height(trialsT), trialsT);
    end
    pseudoSessions_contrast(iP) = IBL_permSession(trialsT,'contrastDiff','pairs');
    pseudoSessions_choice(iP)   = IBL_permSession(trialsT,'choice','pairs');
    pseudoSessions_feedback(iP) = IBL_permSession(trialsT,'feedbackType','pairs');
end

if protoIsBiased
    save(pseudoPath, 'pseudoSessions_bias','pseudoSessions_contrast', ...
        'pseudoSessions_choice','pseudoSessions_feedback');
    pseudoSets.pseudoSessions_bias = pseudoSessions_bias;
else
    save(pseudoPath, 'pseudoSessions_contrast','pseudoSessions_choice','pseudoSessions_feedback');
end
pseudoSets.pseudoSessions_contrast = pseudoSessions_contrast;
pseudoSets.pseudoSessions_choice = pseudoSessions_choice;
pseudoSets.pseudoSessions_feedback = pseudoSessions_feedback;
fprintf('Done!\n');
end

