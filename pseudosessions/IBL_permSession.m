function pseudoS = IBL_permSession(trialsT, cond_to_perm, permMethod)
% IBL_permSession Generate a pseudo-session by permuting one trial variable.
%
% pseudoS = IBL_permSession(trialsT, cond_to_perm)
% pseudoS = IBL_permSession(trialsT, cond_to_perm, permMethod)
%
% Creates a pseudo-session in which the values of a specified trial
% variable are reassigned while preserving the distribution of that
% variable within blocks and within combinations of other task variables.
% Permutations are performed independently within each probabilityLeft
% block to avoid introducing across-block correlations.
%
% INPUTS
%   trialsT
%       Table containing trial-by-trial behavioral data.
%
%   cond_to_perm
%       Name of the variable to permute. Supported values:
%           'choice'
%           'feedbackType'
%           'contrastDiff'
%
%   permMethod (optional)
%       Method used to reassign cond_to_perm within each block and
%       condition-combination.
%
%       'random' (default)
%           Randomly permute cond_to_perm across all trials belonging to
%           the same block and condition-combination.
%
%       'pairs'
%           Construct consecutive pairs of trials within each 
%           condition-combination across the full session (ignoring block 
%           boundaries). Within each pair, cond_to_perm is either left 
%           unchanged or swapped with 50% probability.
%
% OUTPUT
%   pseudoS
%       Structure containing the pseudo-session data with the same fields
%       as the selected condition variables.
%
% NOTES
%   - The variables held fixed depend on cond_to_perm (see switch block
%     below).
%   - The 'pairs' method preserves local trial structure more strongly than
%     full random permutation while still breaking trial-by-trial
%     associations involving cond_to_perm.

if nargin < 3 || isempty(permMethod)
    permMethod = 'random';
end

condFields_names = {'probabilityLeft','contrastDiff', 'choice', 'feedbackType'};
condFields_types = {'double','double','double','double'};

switch cond_to_perm
    case 'choice'
        conds_to_fix = {'probabilityLeft','contrastDiff'};
    case 'feedbackType'
        conds_to_fix = {'probabilityLeft','choice'};
        %conds_to_fix = {'probabilityLeft','contrastDiff'};
    case 'contrastDiff'
        conds_to_fix = {'probabilityLeft','choice'};
end

%find unique conditions to fix
trialsT_conds = trialsT(:,conds_to_fix);
[uniqueCondsT,~,trialCondition] = unique(trialsT_conds,'rows');

%pseudoSession = table('Size',[size(trialsT,1) size(condFields_names,2)],'VariableTypes',condFields_types,'VariableNames',condFields_names);
pseudoSessT = trialsT(:,condFields_names);

nTrials = size(trialsT,1);
blockTransitions = find(abs(diff(trialsT.probabilityLeft))>0.2)';
blockFirstTrials = [1 blockTransitions+1];
blockLastTrials = [blockTransitions nTrials];

switch permMethod

    case 'random'

        % Original behavior: permute within block and condition-combination
        for iBlock = 1:length(blockFirstTrials)

            iTrials_thisBlock = false(nTrials,1);
            iTrials_thisBlock(blockFirstTrials(iBlock):blockLastTrials(iBlock)) = true;

            for i = 1:size(uniqueCondsT,1)

                iCnds = find(trialCondition == i & iTrials_thisBlock);

                pseudoSessT(iCnds,cond_to_perm) = ...
                    pseudoSessT(iCnds(randperm(numel(iCnds))),cond_to_perm);

            end
        end

    case 'pairs'

        % Pairwise swapping across the entire session
        for i = 1:size(uniqueCondsT,1)

            iCnds = find(trialCondition == i);

            nPairs = floor(numel(iCnds)/2);

            if nPairs == 0
                continue
            end

            pairIdx = reshape(iCnds(1:2*nPairs),2,[])';

            for iPair = 1:nPairs

                if rand < 0.5
                    pseudoSessT(pairIdx(iPair,:),cond_to_perm) = ...
                        pseudoSessT(fliplr(pairIdx(iPair,:)),cond_to_perm);
                end

            end
        end

    otherwise

        error('Unknown permMethod: %s',permMethod);

end

pseudoS = table2struct(pseudoSessT,'ToScalar',true);