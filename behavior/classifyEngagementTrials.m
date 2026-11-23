function trialsT = classifyEngagementTrials(trialsT,maxRT,nTrialsDisengaged)
%CLASSIFYENGAGEMENTTRIALS Classify NoGo and disengaged trials.
%
% Inputs:
%   trialsT            : IBL trials table
%   maxRT              : maximum response time (default: 2 s)
%   nTrialsDisengaged  : minimum consecutive NoGo trials defining
%                        a disengaged period (default: 3)
%
% Adds two logical columns:
%   noGoTrials         : response time >= maxRT or NaN
%   disengagedTrials   : runs of >= nTrialsDisengaged consecutive NoGo
%
% If sessionID exists, runs cannot cross session boundaries.

if nargin < 2 || isempty(maxRT)
    maxRT = 2;
end
if nargin < 3 || isempty(nTrialsDisengaged)
    nTrialsDisengaged = 3;
end

responseTime = trialsT.feedback_times - trialsT.goCue_times;
noGoTrials = ~(responseTime(:) < maxRT);

n = height(trialsT);
disengagedTrials = false(n,1);

if ismember('sessionID',trialsT.Properties.VariableNames)
    sessionID = trialsT.sessionID;
else
    sessionID = ones(n,1);
end

for s = unique(sessionID(:))'
    ids = find(sessionID == s);
    runEdges = diff([false; noGoTrials(ids); false]);
    runStarts = find(runEdges == 1);
    runEnds = find(runEdges == -1)-1;

    for iRun = 1:numel(runStarts)
        if runEnds(iRun)-runStarts(iRun)+1 >= nTrialsDisengaged
            disengagedTrials(ids(runStarts(iRun):runEnds(iRun))) = true;
        end
    end
end

trialsT.noGoTrials = noGoTrials;
trialsT.disengagedTrials = disengagedTrials;

end