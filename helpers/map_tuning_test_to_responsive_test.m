function idx = map_tuning_test_to_responsive_test(testName, respNames)
% Map a tuning-stat test name to the corresponding taskResponsiveP column.
%
% Rules:
% - if test contains "stimOn"         -> startsWith signrank_stimOn
% - if test contains "choiceMovement" -> startsWith signrank_choiceMovement
% - if test contains "feedback"       -> startsWith signrank_feedback
% - if test contains "block" or "probabilityLeft" -> startsWith ranksum_baseline_block_probabilityLeft
%
% Returns [] if no match found.

testName = string(testName);
respNames = string(respNames);
idx = [];

if isempty(respNames)
    return
end

if contains(testName, "stimOn", 'IgnoreCase', true) && ~contains(testName, "probabilityLeft", 'IgnoreCase', true)
    idx = find(startsWith(respNames, "signrank_stimOn", 'IgnoreCase', true), 1);

elseif contains(testName, "choiceMovement", 'IgnoreCase', true)
    idx = find(startsWith(respNames, "signrank_choiceMovement", 'IgnoreCase', true), 1);

elseif contains(testName, "feedback", 'IgnoreCase', true)
    idx = find(startsWith(respNames, "signrank_feedback", 'IgnoreCase', true), 1);

elseif contains(testName, "block", 'IgnoreCase', true) || contains(testName, "probabilityLeft", 'IgnoreCase', true)
    idx = find(startsWith(respNames, "ranksum_baseline_block_probabilityLeft", 'IgnoreCase', true), 1);
end
end