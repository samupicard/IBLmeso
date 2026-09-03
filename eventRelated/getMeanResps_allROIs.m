function [meanResps, typeVals, validType] = getMeanResps_allROIs( ...
    resps, trialsT, typeFld, typeVals, condFlds, minTrialsPerCombo)

% getMeanResps_all computes mean responses (PETH) per trial type, for a 3D
% matrix of nROIs x nTrials x nTimepoints
%
%   inputs
%     resps     - nROIs x nTrials x nTimepoints, matrix of time-locked responses
%     trialsT   - nTrials x nCols, trials table with trial types and event latencies
%     typeFld   - string. Trial type to group by (e.g. 'contrastDiff')
%     typeVals  - (optional) values of typeFld to select
%
%   additional input for optional balanced averaging across condition-combinations:
%
%     condFlds: cellstr of additional trial fields (e.g. {'choice','probabilityLeft','feedback'}).
%     minTrialsPerCombo: scalar (default: 10)
%
%   If provided, for each type value we:
%     (1) compute mean within each unique combination of condFlds
%     (2) grand-average across those per-combination means (equal weight per combo),
%         ignoring empty combos.
%
%   outputs
%     avgResps  - nTypeVals x nTimepoints matrix. Mean response for each
%                 trial type specified in typeVals.
%     typeVals  - vector of trial-type values corresponding to rows of
%                 avgResps.
%     validType: 1 x nType logical, true if that type had ALL combos present with >=k trials
%
%   Samuel Picard (Jan 2026)

if nargin < 4, typeVals = []; end
if nargin < 5, condFlds = []; end
if nargin < 6 || isempty(minTrialsPerCombo), minTrialsPerCombo = 2; end

vals = trialsT{:, typeFld};
if isempty(typeVals)
    typeVals = unique(vals);
end

% map trials to main type groups 1..K, keep only selected types
[tfMain, gMain] = ismember(vals, typeVals);
K = numel(typeVals);

resps = resps(:, tfMain, :);
tSel  = trialsT(tfMain, :);
gMain = gMain(tfMain);

nROIs    = size(resps,1);
nTrialsK = size(resps,2);
nT       = size(resps,3);

% reshape responses to trials x (ROIs*time)
X = resps;
validSamp = ~isnan(X);
X(~validSamp) = 0;

X2 = reshape(permute(X,[2 1 3]), nTrialsK, []);
V2 = reshape(permute(validSamp,[2 1 3]), nTrialsK, []);

X2 = double(X2);
V2 = double(V2);

% ---- CASE 1: no balancing requested -> standard trial-weighted mean ----
if isempty(condFlds)
    S = sparse(1:nTrialsK, gMain(:), 1, nTrialsK, K);
    sumByType   = S' * X2;
    countByType = S' * V2;

    A = sumByType ./ countByType;
    A(countByType==0) = NaN;

    meanResps = permute(reshape(A, K, nROIs, nT), [2 1 3]); % nROIs x K x nT
    validType = true(1,K);
    validType(countByType(:,1)<minTrialsPerCombo) = false;
    return
end

% ---- CASE 2: strict balanced mean across condition-combinations ----
% condition-combo id per trial (global across selected trials)
try
    condTbl = tSel(:, condFlds);
catch
    error('One or more fields in condFlds were not found in trialsT.');
end

% -------------------------------------------------------------------------
% Exclude rare conditions from BALANCING ONLY:
%   - choice == 0
%   - feedbackType == 0
%
% These trials are ignored when defining condition-combos and when checking
% combo completeness / balanced averaging.
% -------------------------------------------------------------------------
balanceMask = true(height(tSel),1);

if ismember('choice', tSel.Properties.VariableNames)
    balanceMask = balanceMask & (tSel.choice ~= 0);
end

if ismember('feedbackType', tSel.Properties.VariableNames)
    balanceMask = balanceMask & (tSel.feedbackType ~= 0);
end

% If nothing is left to balance on, return NaNs / invalid
if ~any(balanceMask)
    meanResps = nan(nROIs, K, nT);
    typeVals  = typeVals(:)';
    validType = false(1,K);
    return
end

% Define combos only on balance-eligible trials
condTbl_bal = condTbl(balanceMask, :);
c_bal = findgroups(condTbl_bal);      % 1..C over eligible trials only
C = max(c_bal);

% Full-length combo vector; excluded trials get combo 0
c = zeros(nTrialsK,1);
c(balanceMask) = c_bal;

% Combined group index: (type, combo) -> 1..K*C, only for eligible trials
gComb_bal = (gMain(balanceMask)-1)*C + c_bal;
KC = K*C;

Scomb = sparse(find(balanceMask), gComb_bal(:), 1, nTrialsK, KC);

% --- STRICT REQUIREMENTS: each type must have ALL eligible combos with >=k trials ---
nTrials_tc = full(sum(Scomb, 1));               % 1 x (K*C)
nTrials_tc = reshape(nTrials_tc, [C, K]);       % C x K  (#trials per combo per type)
validType  = all(nTrials_tc >= minTrialsPerCombo, 1);   % 1 x K

if ~all(validType)
    meanResps = nan(nROIs, K, nT);
    typeVals  = typeVals(:)';
    validType = false(1,K);
    return
end

% First-stage: mean within each (type, combo)
sum_tc   = Scomb' * X2;       % (K*C) x P
count_tc = Scomb' * V2;       % (K*C) x P

mean_tc = sum_tc ./ count_tc;
mean_tc(count_tc==0) = NaN;

% Reshape to C x K x P (P = nROIs*nT)
P = size(X2,2);
mean_tc = reshape(mean_tc, [C, K, P]);          % C x K x P

% % Second-stage: equal-weight mean across eligible combos, but ONLY for valid types
% A = nan(K, P);
% if any(validType)
%     A(validType, :) = squeeze(nanmean(mean_tc(:, validType, :), 1));  % (#valid) x P
% end

% Second-stage: equal-weight mean across all combos
A = squeeze(mean(mean_tc, 1, 'omitnan'));       % K x P

% back to nROIs x K x nT
meanResps = permute(reshape(A, K, nROIs, nT), [2 1 3]);

end
